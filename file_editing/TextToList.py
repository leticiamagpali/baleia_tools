def create_list(list_path):
    '''This function takes a text file converts it into a list'''

    list = []
    with open(list_path, "r", encoding="utf-8") as text_file:
        list = text_file.read().splitlines()
    return list
