import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.transforms as tran
from collections import OrderedDict
import sys
#read in obs,mt1_on2, mt1_tim2,mt3,mt5
def get_chi_val(obs_dataset, other_dataset):
	sum=0.0
	for species in obs_dataset.storage:
		for line1, ew1 in obs_dataset.storage[species]:
			for line2, ew2 in other_dataset.storage[species]:
				if line1==line2:
					sum+=float(float((ew1-ew2)**2)/ew2)
					break
	return round(sum,3)		
def get_chi_val_interp(obs_dataset, other_dataset1, other_dataset2):
	sum=0.0
	for species in obs_dataset.storage:
		for line1, ew1 in obs_dataset.storage[species]:
			for line2, ew2 in other_dataset1.storage[species]:
				if line1==line2:
					for line3, ew3, in other_dataset2.storage[species]:
						if line2==line3:
							obs=ew1
							ex=float((ew2+ew3)/2)
							sum+=float(float((obs-ex)**2)/ex)
							break
					break
	return round(sum,3)		

def read_depth_data(filename,lines):
	f = open(filename, "r")
	f_lines = f.readlines()
	f.close()
	storage = dict()
	read_cols=[6,1,2]
	for l in f_lines:
		line = l.split()
		if '(' in line[6]:
			k = line[6][:line[6].index('(')]
			if not k in storage:
				storage[k]=dict()
			if float(line[1]) in lines:
				storage[k][float(line[1])]=float(line[2])
			if not lines:
				storage[k][float(line[1])]=float(line[2])
	return storage
class Dataset:
	def __init__(self, filename, marker, color, label,lines):
		self.storage=self.read_data(filename,lines)
		self.marker=marker
		self.color=color
		self.alt_color=(color[0]/2,color[1]/2,color[2]/2)
		self.label=label

	def read_data(self, filename,lines):
		f = open(filename, "r")
		f_lines = f.readlines()
		f.close()
		storage = dict()
		read_cols=[13,12,4]#12-Line location, 4-EW, 13-Transition
		for l in f_lines[4:]:
			line = l.split()
			if len(line)>13 and '(' in line[13]:
				k = line[13][:line[13].index('(')]
				if not k in storage:
					storage[k]=list()
				if float(line[12]) in lines:
					storage[k].append((float(line[12]),float(line[4])))
				if not lines:
					storage[k].append((float(line[12]),float(line[4])))
		return storage
add_ins=['_SI','_O','_C','_N']
files = ['EW_DATA_O_LINES','EW_DATA_MT1_ON2','EW_DATA_MT1','EW_DATA_MT1_TIM2','EW_DATA_MT3_ON2','EW_DATA_MT3','EW_DATA_MT3_TIM2','EW_DATA_MT5_ON2','EW_DATA_MT5','EW_DATA_MT5_TIM2']
markers=['x','s','s','s','s','s','s','s','s','s']
colors=[(0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)]
labels=[ 'OBS','MT1','MT3','MT5','MT1','MT3','MT5','MT1','MT3','MT5']
count=0
data_storage=list()
data_storage.append(Dataset(files[0],markers[count],colors[count],labels[count],[]))
count+=1
obs_data=data_storage[0]
obs_lines=[]
for species in obs_data.storage:#reads in obs lines
	for line, ew in obs_data.storage[species]:
		obs_lines.append(line)

for f in files[1:]:#reads in all lines into storage
	data_storage.append(Dataset(f,markers[count],colors[count],labels[count],obs_lines))
	count+=1


depth_storage=read_depth_data("LINE_ID",obs_lines)

depth_data_storage=list()
for data in data_storage:
	new_storage=dict()
	for species in data.storage:
		for line, ew in data.storage[species]:
			if not species in new_storage:
				new_storage[species]=list()
			if line in obs_lines:
				new_storage[species].append((depth_storage[species][line],ew))
	depth_data_storage.append(new_storage)

depth_table=dict()

obs_depth_data=depth_data_storage[0]
plt_list=sorted(obs_depth_data["O2"], key=lambda tup: tup[0])
#plt_list2=sorted(obs_depth_data["OIII"], key=lambda tup: tup[0])

for op,ew in plt_list:
	if not op in depth_table:
		depth_table[op]=list()
	depth_table[op].append(ew)

#plt.plot(*zip(*plt_list),color=obs_data.color,label=obs_data.label)
#plt.plot(*zip(*plt_list2),color="grey",label=obs_data.label+'-OIII')

for i in range(1,len(data_storage)):
	plt_list=sorted(depth_data_storage[i]["O2"], key=lambda tup: tup[0])
	for op,ew in plt_list:
		if not op in depth_table:
			depth_table[op]=list()
		depth_table[op].append(ew)
	#plt.plot(*zip(*plt_list),color=data_storage[i].color,label=data_storage[i].label)
	#plt.plot(*zip(*depth_data_storage[i]["OIII"]),color=data_storage[i].color,label=data_storage[i].label)
fig, ax = plt.subplots()	
#plt.xlabel(r'Optical Depth($\tau$)')
#plt.ylabel(r'Equivalent Width(m$\AA$)')
#plt.xlim([0,7.5])
#plt.ylim([0,165])
#handles, labels = plt.gca().get_legend_handles_labels()
#by_label = OrderedDict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys())
# hide axes
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')

txt_table=[[],[],[],[],[]]

for i in range(len(files[1:])):
	chi_val = get_chi_val(data_storage[0],data_storage[i+1])
	txt_table[int(i/3)*2].append(chi_val)
txt_table[1].append(get_chi_val_interp(data_storage[0],data_storage[1],data_storage[4]))
txt_table[1].append(get_chi_val_interp(data_storage[0],data_storage[2],data_storage[5]))
txt_table[1].append(get_chi_val_interp(data_storage[0],data_storage[3],data_storage[6]))
txt_table[3].append(get_chi_val_interp(data_storage[0],data_storage[4],data_storage[7]))
txt_table[3].append(get_chi_val_interp(data_storage[0],data_storage[5],data_storage[8]))
txt_table[3].append(get_chi_val_interp(data_storage[0],data_storage[6],data_storage[9]))
