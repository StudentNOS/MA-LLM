from openai import OpenAI

import os
client = OpenAI()

def generate_pddl(parsed_data):
    """ Generate a patient-directed discharge letter. """
    prompt_data = "\n".join(f"{key}: {value}" for key, value in parsed_data.items())
    
    response = client.chat.completions.create(
    model="gpt-4-1106-preview",
    messages=[
        {"role": "system", "content": PROMPT_LETTER},
        {"role": "user", "content": prompt_data}],
    max_tokens=2048
    )
    return response.choices[0].message.content