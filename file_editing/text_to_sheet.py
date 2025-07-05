import openpyxl


def text_to_spreadsheet(input_file, output_file):
    # Create a new Excel workbook
    workbook = openpyxl.Workbook()
    sheet = workbook.active

    # Open the text file for reading
    with open(input_file, 'r') as file:
        # Read each line from the text file and add it to the Excel sheet
        for row_index, line in enumerate(file, start=1):
            # Assuming each line is a separate string; you may need to adjust this logic
            sheet.cell(row=row_index, column=1, value=line.strip())

    # Save the workbook to an Excel file
    workbook.save(output_file)


# Replace 'input.txt' with the path to your text file and 'output.xlsx' with the desired output Excel file name
text_to_spreadsheet('input.txt', 'output.xlsx')
