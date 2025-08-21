import random

def scramble():
    # Read numbers from the input file
    with open("Initial.txt", "r") as file:
        numbers = file.read().splitlines()
    
    # Shuffle the numbers randomly
    random.shuffle(numbers)
    
    # Write the shuffled numbers to the output file
    with open("Initial.txt", "w") as file:
        file.write("\n".join(numbers))

    print("Numbers have been shuffled and saved.")