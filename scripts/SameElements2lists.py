# !/ usr/bin/ python3

'''
Author Julia Lienard
Date 2024-02-28

Description: Python program to find the common elements in two different lists
Usage: python SameElements2lists.py file1 file2 output
'''

import sys

List1 = sys.argv[1]
List2 = sys.argv[2]
CommonList = sys.argv[3]

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    return (a_set & b_set)

List1ID = []
List2ID = []

with open(List1, "r") as file1, open(List2) as file2, open(CommonList, "w") as output:
    for line in file1:
        List1ID.append(line.strip())
    for line in file2:
        List2ID.append(line.strip())
    common_element = common_member(List1ID, List2ID)
    for element in common_element:
        output.write(f'{element}\n')
