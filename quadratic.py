#this program solves a quadratic equation ax^2+bx+c=0
class Quadratic:
    def __init__(self,a,b,c):
        self.a=a
        self.b=b
        self.c=c
        print("The quadratic equation is:", a,"x^2 + ",b,"x + ",c)
        
    def solutions(self):
        discriminant=self.b*self.b-4*self.a
        print("discriminant = ",discriminant)
        if (discriminant>0):
            root1=(-self.b-discriminant**0.5)/(2*self.a)
            root2=(-self.b+discriminant**0.5)/(2*self.a)
            print(root1)
            print(root2)
            roots=[root1,root2]
            return roots
        elif (discriminant==0):
            root=-self.b/(2*self.a)
            print(root)
            roots=[root]
            return roots
        else:
            print("Since the discriminant is negative, there are no real solutions.")
            root1=(-self.b-discriminant**0.5)/(2*self.a)
            root2=(-self.b+discriminant**0.5)/(2*self.a)
            #print(root1)
            #print(root2)
            roots=[root1,root2]
            return roots
                    
import numpy as np
a=np.random.randint(1,10)
b=np.random.randint(1,10)
c=np.random.randint(1,10)
qd2=Quadratic(a,b,c)
qd2.solutions()
