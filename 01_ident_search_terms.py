# Input: Study protocol
# Output: search terms

import openai

def generate_openai_response(prompt):
    openai.api_key = 'API-key upon request'

    response = openai.Completion.create(
      engine="gpt-3.5-turbo",  # You can choose different models
      prompt=prompt,
      max_tokens=150  # You can adjust the number of tokens
    )

    return response.choices[0].text.strip()

prompt = "Analyze the following study plan: study plan"
search_tearms = "From this study plan, generate a list of possible search terms..."
format = "Specify format of response"
response = generate_openai_response(prompt, search_terms, format)
print(response)
