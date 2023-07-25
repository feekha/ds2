# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:53:49 2022

@author: HENROG
"""

# class Person:
#   # def __init__(self, fname, lname):
#   #   self.firstname = fname
#   #   self.lastname = lname

#   def printname(self, fname, lname):
#     self.val = 200  
#     print(fname, lname)
#     print(self.graduationyear)

# class Student(Person):
#   def __init__(self, fname, lname, year):
#     # super().printname(fname, lname)
#     self.graduationyear = year
    
#   def do(self, gradyear):
#       print(gradyear)
      

# x = Student("Mike", "Olsen", 2019)
# # print(x.graduationyear)
# x.printname(fname="Hendrik", lname="Rogoll")


import pandas as pd


df = pd.DataFrame({'A': [1, 2, 3, 4, 5],
                   'B': [6, 7.75, 8, 9, 10],
                   'C': [11, 12, 13, 14, 15],
                   'D': [16, 17, 18, 19, 20]},
                  index=pd.date_range('2022-01-01', periods=5, freq='D'))

# Create bins
bins = pd.interval_range(start=3.75, end=16.25, freq=0.5, closed='right')

# Use cut function to create the bins
binned_dfs = []
for col in df.columns:
    df[str(col)+'_binned'] = pd.cut(df[col], bins=bins) 
    df_binned = df[[col,str(col)+'_binned']].groupby(str(col)+'_binned').count()
    binned_dfs.append(df_binned)

df_binned = pd.concat(binned_dfs, axis = 1)

# # Create a new categorical index with the same categories as the original index
# new_index = pd.CategoricalIndex(categories, categories)

# # Set the new categorical index as the index of the sorted dataframe
# new_df = df_grouped.set_index(new_index)

# # Remove the 'bins' column
# new_df = new_df.drop(columns=['bins'])  