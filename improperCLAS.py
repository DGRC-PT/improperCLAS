from sys import argv
import statistics

def is_reverse(v):
	if int(v) & 16:
		return True
	else:
		return False

def is_proper(v):
	if int(v) & 2:
		return True
	else:
		return False

def estimate_isize(infile):
	f=open(infile)
	m=[]
	for el in f:
		if el.startswith("@")==False:
			line=el.split("\t")
			if is_proper(line[1])==True and float(line[8])>0:
				m.append(float(line[8]))
			if len(m)==10000:
				break
	f.close()
	print("IS median: "+str(statistics.median(m)))
	print("IS stdev: "+str(statistics.stdev(m)))
	return statistics.median(m), statistics.stdev(m)

def is_unmaped(v):
	if int(v) & 4 or int(v) & 8:
		return True
	else:
		return False

def read_sam_for_CNVs(infile, me, st):#infile - read sorted sam file
	if me=="" and st=="":
		media, std=estimate_isize(infile)
	else:
		print("inputed IS median: "+me)
		print("inputed stdev: "+st)
		media=float(me)
		std=float(st)
	lim=media+(3*std)
	f=open(infile)
	name=""
	dell=[]
	dup=[]
	inv=[]
	trans=[]
	for i in f:
		if i.startswith("@")==False:
			line=i.split("\t")
			if name=="" and is_unmaped(line[1])==False:
				name=line[0]
				temp=line
			elif line[0]==name and is_unmaped(line[1])==False:
				if line[2]==temp[2]:
					if int(temp[3])<int(line[3]) and is_reverse(temp[1])==True and is_reverse(line[1])==False:
						if abs(int(temp[8]))>int(lim):
							dup.append(name+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t26\t"+temp[4]+"\t26\t"+line[4]+"\t"+temp[9]+"\t"+temp[10].strip()+"\t"+line[9]+"\t"+line[10].strip())
					if int(temp[3])>int(line[3]) and is_reverse(temp[1])==False and is_reverse(line[1])==True:
						if abs(int(temp[8]))>int(lim):
							dup.append(name+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t26\t"+line[4]+"\t26\t"+temp[4]+"\t"+line[9]+"\t"+line[10].strip()+"\t"+temp[9]+"\t"+temp[10].strip())
					if int(temp[3])<int(line[3]) and is_reverse(temp[1])==False and is_reverse(line[1])==True:
						if abs(int(temp[8]))>int(lim):
							dell.append(name+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t26\t"+temp[4]+"\t26\t"+line[4]+"\t"+temp[9]+"\t"+temp[10].strip()+"\t"+line[9]+"\t"+line[10].strip())
					if int(temp[3])>int(line[3]) and is_reverse(temp[1])==True and is_reverse(line[1])==False:
						if abs(int(temp[8]))>int(lim):
							dell.append(name+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t26\t"+line[4]+"\t26\t"+temp[4]+"\t"+line[9]+"\t"+line[10].strip()+"\t"+temp[9]+"\t"+temp[10].strip())
				name=""
				temp=[]
			else:
				name=""
				temp=[]
	f.close()
	return dell, dup

def sam_for_other(infile):
	ch=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y", "MT"]
	f=open(infile)
	name=""
	inv=[]
	trans=[]
	for i in f:
		if i.startswith("@")==False:
			line=i.split("\t")
			if len(line[2])>3:
				line[2]=line[2][3:]
				line[6]=line[6][3:]
			if name=="" and is_unmaped(line[1])==False:
				name=line[0]
				temp=line
			elif line[0]==name and is_unmaped(line[1])==False:
				if line[2]!=temp[2]:
					if ch.index(temp[2])<ch.index(line[2]):
						trans.append(name+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t26\t"+temp[4]+"\t26\t"+line[4]+"\t"+temp[9]+"\t"+temp[10].strip()+"\t"+line[9]+"\t"+line[10].strip())
					if ch.index(temp[2])>ch.index(line[2]):
						trans.append(name+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t26\t"+line[4]+"\t26\t"+temp[4]+"\t"+line[9]+"\t"+line[10].strip()+"\t"+temp[9]+"\t"+temp[10].strip())
				elif is_reverse(line[1])==is_reverse(temp[1]):
					if int(temp[3])<int(line[3]):
						inv.append(name+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t26\t"+temp[4]+"\t26\t"+line[4]+"\t"+temp[9]+"\t"+temp[10].strip()+"\t"+line[9]+"\t"+line[10].strip())
					else:
						inv.append(name+"\t"+line[2]+"\t"+line[3]+"\t"+line[1]+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[1]+"\t26\t"+line[4]+"\t26\t"+temp[4]+"\t"+line[9]+"\t"+line[10].strip()+"\t"+temp[9]+"\t"+temp[10].strip())
				name=""
				temp=[]
			else:
				name=""
				temp=[]
	f.close()
	return inv, trans


def write(l, outfile):
	out=open(outfile, "w")
	for el in l:
		out.write(el+"\n")
	out.close()

if len(argv)>6:
	m=argv[6]
	st=argv[7]
else:
	m=""
	st=""

#python subbamstat.py samfile, del_out, dup_out, inv_out, trans_out
dell, dup=read_sam_for_CNVs(argv[1], m, st)
write(dell, argv[2])
write(dup, argv[3])
dell=[]
dup=[]
inv, trans=sam_for_other(argv[1])
write(inv, argv[4])
write(trans, argv[5])

