import sys
with open(sys.argv[1],"r",encoding="utf-8") as f:
    lines = f.readlines()
    #print(lines)
with open(sys.argv[2],"w",encoding="utf-8") as f_w:
	i = 0
	for line in lines:
		if '# Name "object23"' in line or '# Name "relief"' in line or '# Name "holes"' in line or '# Name "outside01"' in line or '# Name "parapet"' in line or '# Name "round_hole"' in line or '# Name "walls"' in line or '# Name "windows"' in line or '# Name "doors"' in line:
			i += 1
			continue
		if i != 0:
			i += 1 
			if i == 9:
				i = 0
				continue
			continue
		f_w.write(line)


