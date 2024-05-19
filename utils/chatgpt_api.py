from openai import OpenAI

def chat(prompt, api_key, model="gpt-3.5-turbo"):
    client = OpenAI(api_key=api_key)

    response = client.chat.completions.create(
        model="gpt-3.5-turbo",  
        messages=[{"role": "system", "content": prompt}],
        max_tokens=2048  
    )

    return response.choices[0].message.content


#response = chat('Say goodbye', 'api-key on request')
#print(response)