# !/ usr/bin/ python3

'''
2024-02-28
refined by ChatGPT, based on my existig SameElements2lists.py script

Python program to find the common elements in different lists

'''
import sys

def common_member(*lists):
    common_set = set(lists[0])
    for lst in lists[1:]:
        common_set &= set(lst)
    return common_set

# Check if at least two list files are provided
if len(sys.argv) < 4:
    print("Usage: python script.py List1.txt List2.txt CommonList.txt [List3.txt ...]")
    sys.exit(1)

list_files = sys.argv[1:]
output_file = list_files.pop()  # Pop the last file for the output

list_data = []

for list_file in list_files:
    with open(list_file, "r") as file:
        current_list = [line.strip() for line in file]
        list_data.append(current_list)

common_elements = common_member(*list_data)

with open(output_file, "w") as output:
    for element in common_elements:
        output.write(f'{element}\n')
