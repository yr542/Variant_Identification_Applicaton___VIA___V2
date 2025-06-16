# Family class. Stores the strings and phenotype for every member of the family.

class Family:
    def __init__(self, ID):
        self.ID = ID
        self.people = []
        self.siblings = []
        self.mother = Person("","","")
        self.father = Person("","","")
        self.child = Person("","","")
        self.hasFather = False
        self.hasMother = False
        self.HPO = []
        self.genes = {}

# Sibling class. Stores the string, sex, and phenotype for a sibling.
class Person:
    def __init__(self, ID, sex, phen):
        self.ID = ID
        self.sex = sex
        self.phen = phen
        self.affected = self.phen == "Affected"
        self.unaffected = self.phen == "Unaffected"
        self.male = self.sex == "Male"
        self.female = self.sex == "Female"
