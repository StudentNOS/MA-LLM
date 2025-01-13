# %load -r 3:25 _init.py
import pathlib
import sammo
from sammo.runners import OpenAIChat
from sammo.base import Template, EvaluationScore
from sammo.components import Output, GenerateText, ForEach, Union
from sammo.extractors import ExtractRegex
from sammo.data import DataTable
from sammo.mutators import BagOfMutators, InduceInstructions, Paraphrase
from sammo.search import BeamSearch
from sammo.instructions import MetaPrompt, Section, Paragraph, InputData
from sammo.dataformatters import PlainFormatter
from sammo.search_op import one_of
import pandas as pd
import json
import requests
import os
from sklearn.metrics import f1_score

if not "OPENAI_API_KEY_NOS" in os.environ:
    raise ValueError("Please set the environment variable OPENAI_API_KEY'.")

_ = sammo.setup_logger("WARNING")  # we're only interested in warnings for now

runner = OpenAIChat(
    model_id="gpt-3.5-turbo-16k",
    api_config={"api_key": os.getenv("OPENAI_API_KEY_NOS")},
    cache=os.getenv("CACHE_FILE", "cache.tsv"),
    timeout=30,
)

# %load -s load_data,accuracy _init.py
def load_data(file_path="./data/train_data.csv"):

    df = pd.read_csv(file_path)
    examples = df.to_dict(orient="records")

    return DataTable.from_records(
        examples,
        input_fields=["Number"],
        output_fields=["Label"],
        constants={"instructions": "Please classify the following examples:"},
    )

mydata = load_data(file_path="./data/train_data.csv")
d_train = mydata.sample(10, seed=42)


def f1(y_true: DataTable, y_pred: DataTable) -> EvaluationScore:
    y_true = y_true.outputs.normalized_values()
    y_pred = y_pred.outputs.normalized_values()
    score = f1_score(y_true, y_pred, average='weighted')

    return EvaluationScore(score)


class InititialCandidates:
    def __init__(self, dtrain):
        self.dtrain = dtrain

    def __call__(self):
        example_formatter = PlainFormatter(
            all_labels=self.dtrain.outputs.unique(), orient="item"
        )

        labels = self.dtrain.outputs.unique()
        instructions = MetaPrompt(
            [
            Paragraph("Instructions: "),
            Paragraph(
                one_of(
                [
                    "Screen the titles of these papers as a good scientist would.",
                    "Screen the titles of these papers as an expert in efficacy and safety of mesenchymal stem cell therapy in liver cirrhosis would.",
                    "Screen these papers to ensure only the most relevant studies on efficacy and safety of mesenchymal stem cell therapy in liver cirrhosis are included.",
                    "Screen these papers to make sure no studies on efficacy and safety of mesenchymal stem cell therapy in liver cirrhosis are missed.",
                    "Do not exclude any papers unless you are absolutely certain they are irrelevant to the meta-analysis.",
                    "Include all papers unless you are absolutely certain they are irrelevant to the meta-analysis.",
                    "Include all papers unless you are ***absolutely certain*** they are irrelevant to the meta-analysis.",
                    "Screen these papers making sure to include all papers that might be relevant. Be generous with inclusion.",
                    "Select the most fitting papers according to my request.",
                    "Which of these papers should I include if I want to produce the most comprehensive meta-analysis possible?",
                    "Select the titles that suggest that the paper investigates the efficacy and safety of mesenchymal stem cells in the treatment of liver cirrhosis.",
                    "Screen these papers like a human would.",
                ]
                ),
                reference_id="instructions",
            ),
            Paragraph("\n"),
            Paragraph(
                f"Output labels: {', '.join(labels)}\n" if len(labels) <= 10 else ""
            ),
            Paragraph(InputData()),
            Paragraph("Output: "),
            ],
            render_as="raw",
            data_formatter=example_formatter,
        )

        return Output(
            instructions.with_extractor("raise"),
            minibatch_size=1,
            on_error="empty_result",
        )

mutation_operators = BagOfMutators(
    InititialCandidates(d_train),
    InduceInstructions("#instructions", d_train),
    Paraphrase("#instructions"),
    sample_for_init_candidates=False,
)


prompt_optimizer = BeamSearch(
            runner,
            mutation_operators,
            f1,
            maximize=True,
            depth=3,
            mutations_per_beam=2,
            n_initial_candidates=4,
            beam_width=4,
            add_previous=True,
    )
prompt_optimizer.fit(d_train)
prompt_optimizer.show_report()

print(prompt_optimizer.best_prompt)

