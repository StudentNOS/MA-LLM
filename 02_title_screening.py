# Input: Search terms
# Output: Pre-selected papers based on titles

from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='API Key upon request')

# Function to fetch all titles from the database
def get_all_titles():
    """
    Fetches all titles from the database.
    """
    query = "SELECT title FROM study"
    titles = execute_query(query, ENSURE)
    return [title[0] for title in titles]

# Function to create a prompt for OpenAI
def generate_prompt(search_term, titles):
    """
    Generates a prompt for OpenAI based on a search term and a list of titles.
    """
    prompt = f"Search Term: {search_term}\n\nList of titles:\n"
    prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))
    prompt += "\n\nWhich of these titles are relevant to the search term?"
    return prompt

# Function to screen titles using OpenAI API
def screen_titles_with_openai(prompt):
    """
    Screens titles using the OpenAI API.
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
    Main function to screen paper titles using GPT.
    """
    search_term = input("Enter search term: ")
    titles = get_all_titles()
    prompt = generate_prompt(search_term, titles)
    screened_titles = screen_titles_with_openai(prompt)

    print("Screened titles:")
    print(screened_titles)

if __name__ == "__main__":
    main()
