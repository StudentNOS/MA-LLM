# Input: Search terms
# Output: Pre-selected papers based on abstracts

from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='API Key upon request')

# Function to fetch all abstracts from the database
def get_all_abstracts():
    """
    Fetches all abstracts from the database.
    """
    query = "SELECT abstract FROM study"
    abstracts = execute_query(query, ENSURE)
    return [abstract[0] for abstract in abstracts]

# Function to create a prompt for OpenAI
def generate_prompt(search_term, abstracts):
    """
    Generates a prompt for OpenAI based on a search term and a list of abstracts.
    """
    prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
    prompt += "\n".join(f"{i}. {abstract}" for i, abstract in enumerate(abstracts, 1))
    prompt += "\n\nWhich of these abstracts are relevant to the search term?"
    return prompt

# Function to screen abstracts using OpenAI API
def screen_abstracts_with_openai(prompt):
    """
    Screens abstracts using the OpenAI API.
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return ""

# Main function for GPT screening
def main():
    """
    Main function to screen paper abstracts using GPT.
    """
    search_term = input("Enter search term: ")
    abstracts = get_all_abstracts()
    prompt = generate_prompt(search_term, abstracts)
    screened_abstracts = screen_abstracts_with_openai(prompt)

    print("Screened abstracts:")
    print(screened_abstracts)

if __name__ == "__main__":
    main()
