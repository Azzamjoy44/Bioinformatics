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

print(find_alphabet("A*AAB"))
        