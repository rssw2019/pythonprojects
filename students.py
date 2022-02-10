class Students:
    def __init__(self,firstname,lastname,semester):
        self.firstname=firstname
        self.lastname=lastname
        self.semester=semester

student1=Students("Harry", "Trueman", "Fall 2021")



print(student1.firstname+" " +student1.lastname)
