# Input: Study protocol
# Output: Inclusion/Exclusion criteria

from openai import OpenAI

client = OpenAI(api_key=api_key_upon_request)

def generate_openai_response(prompt):

    response = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[{"role": "system", "content": prompt}],
        max_tokens=2048
    )
    
    response = response.choices[0].message.content
    return response

# Example usage
prompt = "Analyze the following study plan: A meta-analysis will be written about the effects of Covid in Asthma patients"
prompt += "From this study plan, generate a list of possible inclusion and exclusion criteria"
prompt += "Specify format of response"
response = generate_openai_response(prompt)
print(response)
