import os

class sample(object):
	def __init__(self, name, total, discordant, concordant, leftover): 
		self.name = name
		self.total = int(total)
		self.unconcordant = int(concordant)
		self.discordant = int(discordant)
		self.concordant = int(concordant)
		self.leftover = int(leftover)

	def __str__(self):
		return str(self.name)
		

def data_clean(data_dic):
 	clean_data = {}
 	flag = False
 	for key,value in data_dic.items():
 		total = value[0].split()[0]
 		unconcordant = value[1].split()[0]
 		concordant_one = value[2].split()[0]
 		condordant_more = value[3].split()[0]
 		concordant = int(condordant_more) + int(concordant_one)
 		discordant = value[6].split()[0]
 		leftover = value[8].split()[0]

 		print("Total",total)
 		print("unconcordant", unconcordant)
 		print("concordant_one", concordant_one)
 		print("condordant_more",condordant_more)
 		print("concordant", concordant) 
 		if (int(total) - int(concordant)) == int(unconcordant):
 			print(key, "passes santiy")
 		else:
 			print("WARNING", key, "does not pass sanity. ")
 			print("Total", total,"minus concordant", concordant,"is not equal to unconcordant", unconcordant)
 			flag = True
 		print("discordant", discordant)
 		print("leftover", leftover)
 		if (int(leftover) + int(discordant)) == int(unconcordant):
 			print(key, "passes sanity 2")
 		else: 
 			print("WARNING", key, "does not pass sanity 2.")
 			flag = True 
 		print("\n")
 		clean_data[key] = (total, unconcordant, concordant, discordant, leftover)

 	if flag == True:
 		print("WARNING one or more samples have not passed the sanity checks. Please look through output.")
 	else:
 		return clean_data

def data_print(data_dic):
	for key,value in data_dic.items(): 
		print("Sample:", key)
		print("Data:", value)
		print("\n")
	return 

def main():

	path = 'C:\\Users\\brendan\\Documents\\sravandevanathan\\run_july_27_2019_gencode\\slurm_err'
	files = []
	sample_data = {}

	# r=root, d=directories, f = files
	for r, d, f in os.walk(path):
		for file in f:
			files.append(os.path.join(r, file))

	#sanity 
	for f in files:
		curr = open(f, "r")
		mark = False
		for line in curr:

			if line[0].isnumeric() == True:
				name = f[80:]
				if name[0] == "2":
					name = name[6:19]
				else: 
					name = name[21:30]
				
				data = curr.read()
				data = data.split("\n")
				sample_data[name] = data
				mark = True

			if mark == True:
				break 
				
		curr.close()

	clean_data = data_clean(sample_data)
	samples = []
	for key,value in clean_data.items():
		samples.append(sample(key, value[0], value[1], value[2], value[3]))

	#ignore unless wants concordant values
	concordant_values = {}
	f = open("concordant.txt", "w")
	for item in samples:
		concordant_values[item.name] = item.concordant
		f.write((item.name +","+ str(item.concordant)+","+ str(item.discordant)+ "," + str(item.leftover)+"\n"))
	print(concordant_values)
main()