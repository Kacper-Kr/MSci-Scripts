import random
umales = 100
imales = 0
ufemales = 99
ifemales = 1
infections = []

def printpopulation():
    print("Infected Males: " +str((imales/(imales+umales))*100)+"%")
    print("Infected Females: " +str((ifemales/(ifemales+ufemales))*100)+"%")

    return ((imales + ifemales)/(imales + ifemales + umales + ufemales))*100
    
def newgeneration(oldumales,oldimales,oldufemales,oldifemales):
    newumales = 0
    newimales = 0
    newufemales = 0
    newifemales = 0
    for i in range(0,ifemales):
        newimales = newimales +1
        newifemales = newifemales +1
    for i in range(0,ufemales):
        newumales =newumales +1
        newufemales =newufemales +1
    newumales = newumales - (1*imales)
    newufemales = newufemales - (1*imales)
    return newumales,newimales,newufemales,newifemales

while ufemales > 0 and len(infections) < 200:
    infected = printpopulation()
    infections.append(infected)
    umales,imales,ufemales,ifemales = newgeneration(umales,imales,ufemales,ifemales)

infected = printpopulation()
infections.append(infected)
print(infections)

def graph():
    y100 =  "100-|"
    y90 =   "    |"
    y80 =   " 80-|"
    y70 =   "    |"
    y60 =   " 60-|"
    y50 =   "    |"
    y40 =   " 40-|"
    y30 =   "    |"
    y20 =   " 20-|"
    y10 =   "    |"
    y0 =    "  0-|"
    for i in range(0,len(infections)):
        y100 += " "
        y90 +=  " "
        y80 +=  " "
        y70 +=  " "
        y60 +=  " "
        y50 +=  " "
        y40 +=  " "
        y30 +=  " "
        y20 +=  " "
        y10 +=  " "
        y0 +=   " "
        if infections[i] >= 10:
            y0 = y0[0:-1] +"#"
        elif infections[i] >= 5:
            y0 = y0[0:-1] +"-"
        elif infections[i] > 0:
            y0 = y0[0:-1] +"_"
        if infections[i] >= 20:
            y10 = y10[0:-1] +"#"
        elif infections[i] >= 15:
            y10 = y10[0:-1] +"-"
        elif infections[i] >= 10:
            y10 = y10[0:-1] +"_"
        if infections[i] >= 30:
            y20 = y20[0:-1] +"#"
        elif infections[i] >= 25:
            y20 = y20[0:-1] +"-"
        elif infections[i] >= 20:
            y20 = y20[0:-1] +"_"
        if infections[i] >= 40:
            y30 = y30[0:-1] +"#"
        elif infections[i] >= 35:
            y30 = y30[0:-1] +"-"
        elif infections[i] >= 30:
            y30 = y30[0:-1] +"_"
        if infections[i] >= 50:
            y40 = y40[0:-1] +"#"
        elif infections[i] >= 45:
            y40 = y40[0:-1] +"-"
        elif infections[i] >= 40:
            y40 = y40[0:-1] +"_"
        if infections[i] >= 60:
            y50 = y50[0:-1] +"#"
        elif infections[i] >= 55:
            y50 = y50[0:-1] +"-"
        elif infections[i] >= 50:
            y50 = y50[0:-1] +"_"
        if infections[i] >= 70:
            y60 = y60[0:-1] +"#"
        elif infections[i] >= 65:
            y60 = y60[0:-1] +"-"
        elif infections[i] >= 60:
            y60 = y60[0:-1] +"_"
        if infections[i] >= 80:
            y70 = y70[0:-1] +"#"
        elif infections[i] >= 75:
            y70 = y70[0:-1] +"-"
        elif infections[i] >= 70:
            y70 = y70[0:-1] +"_"
        if infections[i] >= 90:
            y80 = y80[0:-1] +"#"
        elif infections[i] >= 85:
            y80 = y80[0:-1] +"-"
        elif infections[i] >= 80:
            y80 = y80[0:-1] +"_"
        if infections[i] >= 100:
            y90 = y90[0:-1] +"#"
        elif infections[i] >= 95:
            y90 = y90[0:-1] +"-"
        elif infections[i] >= 90:
            y90 = y90[0:-1] +"_"
        if infections[i] >= 110:
            y100 = y100[0:-1] +"#"
        elif infections[i] >= 105:
            y100 = y100[0:-1] +"-"
        elif infections[i] >= 100:
            y100 = y100[0:-1] +"_"
        
    print(y100 +"\n" +y90 +"\n"  +y80 +"\n"  +y70 +"\n"  +y60 +"\n"  +y50 +"\n"  +y40 +"\n"  +y30 +"\n"  +y20 +"\n"  +y10 +"\n"  +y0)
print("")
graph()
print("")
print("Population reached " +str(infections[-1])+"% infection after " +str(len(infections)-1)+" generations.")
        
