crs = open("Iter_1.txt", "r")
iter = 0;
for columns in ( raw.strip().split() for raw in crs ):  
    iter = iter+1;
matrix = [[0]*3 for i in range(iter)]
crs.closed
crs = open("Iter_1.txt", "r")
iter = 0;
for columns in ( raw.strip().split() for raw in crs ):
    matrix[iter][0]=int(columns[0])
    matrix[iter][1]=int(columns[1])
    matrix[iter][2]=float(columns[2])
    iter = iter+1;
crs.closed


