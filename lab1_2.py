s = "ACGGGCATATGCGC"

alphabet = []

def find_alphabet(string):
    for char in string:
        char_exists = False
        for symbol in alphabet:
            if char == symbol:
                char_exists = True
                break
        if not char_exists:
            alphabet.append(char)
    return alphabet

print(find_alphabet(s))

def find_percentages(string):
    percentages = {}
    length = len(string)
    for char in alphabet:
        count = 0
        for symbol in string:
            if char == symbol:
                count += 1
        percentages[char] = (count / length) * 100
    return percentages

print(find_percentages(s))
