**1.	Join GitHub repository for access to all relevant scripts** 
(https://github.com/BerlinExchangeMedicine/ensure)


**2.	Load utils script “dbconnect.py“**

-> creates SQLite-database “ensure.sqlite” on your PC 



**3.	Load utils script “PubMed_for_ID_Titles_Abstracts.py**

-> asks search terms for PubMed query

-> searches PubMed API for the search terms you enter

-> saves IDs to SQLite-database

-> creates Excel-table from database (“ensure.xlsx”)

-> outputs how many papers are in database
![image](https://github.com/BerlinExchangeMedicine/ensure/assets/133876003/c3e68c1e-d0ed-4e4a-9f3e-eff9cb736656)


**4.	 Load script “title_screening.py”**

-> requires OpenAI API-key of NOS-account (can be given upon request)

-> requires you to download the “initial_search.txt” and “golden_ref.txt” files (soon on GitHub)

-> contains prompt for prompt-engineering

-> asks search terms for screening

-> screens the extracted PubMed titles for the search terms you enter 
![image](https://github.com/BerlinExchangeMedicine/ensure/assets/133876003/e31850a7-8c9d-4f3d-a8a7-4cdf06943d74)

_-> Outputs 3 percentages:_ 

a.	**Precision:** GPT x Gold Standard (i.e., What percentage of the IDs in the gold standard did GPT also find?) = should be as close to 100% as possible 

b.	**Proportion covered (Gold standard):** Gold Standard X Initial Search (i.e., What percentage of the initial search IDs was extracted manually in gold standard) = serves as reference value for c

c.	**Proportion covered (GPT):**GPT X Initial Search (What percentage of the initial search IDs was extracted by GPT) = should be as close to reference value/b as possible 
![image](https://github.com/BerlinExchangeMedicine/ensure/assets/133876003/a38c3505-d596-450f-85de-fd6b2e7a9bdf)


**5.	Load script “abstract_screening.py”**

-> Ideally outputs the same as “title_screening.py”

-> If not, gives GPT-output for manual comparison of extracted abstracts 

-> Only change content prompt and leave format prompts unchanged
![image](https://github.com/BerlinExchangeMedicine/ensure/assets/133876003/697aa131-4096-4329-8373-63dc8b69a4db)
