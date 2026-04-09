import csv 
import sys 
import argparse 


def compare_csv(file1, file2):
    with open(file1, 'r', encoding='utf-8') as f1, open(file2, 'r', encoding='utf-8') as f2:
        for r, (a, b) in enumerate(zip(csv.reader(f1), csv.reader(f2))):
            for c, (x, y) in enumerate(zip(a, b)):
                if x != y:
                    print(f'Row {r+1}, Col {c+1}: {x!r} -> {y!r}')


def main():
    parser = argparse.ArgumentParser(description='Compare two CSV files.')
    parser.add_argument('file1', help='The first CSV file.')
    parser.add_argument('file2', help='The second CSV file.')
    args = parser.parse_args()
    compare_csv(args.file1, args.file2)


if __name__ == '__main__':
    main()