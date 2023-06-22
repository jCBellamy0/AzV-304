import matplotlib.pyplot as plt
import sys
def read_data(filename):
	f = open(filename, "r")
	f_lines = f.readlines()
	f.close()
	storage = dict()
	read_cols=[13,2,4]#2-Line location, 4-EW, 13-Transition
	for l in f_lines[4:]:
		line = l.split()
		if '(' in line[13]:
			k = line[13][:line[13].index('(')]
			if not k in storage:
				storage[k]=list()
			storage[k].append((float(line[2]),float(line[4])))
	return storage
#filename = sys.argv[1]
files = sys.argv[1:]
mt_vals=[1,3,5]
mt_storage=list()
for f_i in range(len(files)):
	temp_str=read_data(files[f_i])
	mt_storage.append(temp_str)
for store in mt_storage[1:]:#normalize to mt=1
	for i in range(len(store["SkIII"])):
		store["SkIII"][i]=(store["SkIII"][i][0],store["SkIII"][i][1]/mt_storage[0]["SkIII"][i][1])
	for i in range(len(store["SkIV"])):
		store["SkIV"][i]=(store["SkIV"][i][0],store["SkIV"][i][1]/mt_storage[0]["SkIV"][i][1])

for i in range(len(mt_storage[0]["SkIII"])):
	mt_storage[0]["SkIII"][i]=(mt_storage[0]["SkIII"][i][0],1.0)
for i in range(len(mt_storage[0]["SkIV"])):
	mt_storage[0]["SkIV"][i]=(mt_storage[0]["SkIV"][i][0],1.0)

mt_ct=0
color_array=['b','g','r','c','y','m']
col_ct=0
for store in mt_storage:
	plt.scatter(*zip(*store["SkIII"]),marker='o',color=color_array[col_ct],label="SkIII-MT="+str(mt_vals[mt_ct]))
	mt_ct+=1
	col_ct+=1
mt_ct=0
for store in mt_storage:
	plt.scatter(*zip(*store["SkIV"]),marker='o',color=color_array[col_ct],label="SkIV-MT="+str(mt_vals[mt_ct]))
	mt_ct+=1
	col_ct+=1


#plt.plot(*zip(*storage["O2"]),label="O2")
plt.xlabel("Center Wavelength (Angstroms)")
plt.ylabel("EW (mA)")
plt.legend()
plt.show()
	