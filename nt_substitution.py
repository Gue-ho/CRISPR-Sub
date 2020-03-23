import time, operator, argparse, xlsxwriter,  sys
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.stats import norm
import numpy as np
from math import pi, log


def rev_comp(seq): return seq.translate(seq.maketrans('ATGC','TACG'))[::-1]

def func(x,a,b,c,d): return a+c*np.exp(-pow((x-b),2)/(2*(pow(d,2))))

class NtSubstitution:

	def __init__(self, args):

		self.wt_seq = args.wt_seq.upper()
		self.t_seq = args.target_seq.upper()
		self.input_file = {'data': args.data_file, 'con': args.con_file}
		self.out_file = args.out_file
		#test option
		self.end_range = 60
		self.marker = 5
		self.filt_r = 5
		self.start_pos = -1
		self.end_pos = -1
		self.line_dict = {'data': {}, 'con': {}}
		self.cnt_dict = {'data':{'all': 0 ,'indi': 0, 'len': 0, 'ins': 0, 'del': 0, 'sub': 0}, 'con':{'all': 0 ,'indi': 0, 'len': 0, 'ins': 0, 'del': 0}}
		self.rate_dict = {'data': {}, 'con': {}}
		self.value_dict = {'data': {}, 'con': {}}

	def seq_env(self):

		print('Set sequence options')

		s_num = 0
		e_num = 0

		for i in range(22, len(self.wt_seq) - 23):

			if self.wt_seq[i:i+23] == self.t_seq:
				self.start_pos = i + 17 - self.end_range
				self.end_pos = i + 17 +self.end_range

				if self.start_pos < 0:
					s_num = (-1) * self.start_pos
					start_pos = 0

				if self.end_pos > len(self.wt_seq):
					e_num = self.end_pos - len(self.wt_seq)
					end_pos = len(self.wt_seq)

				self.s_seq = self.wt_seq[i + 17 - self.filt_r: i + 17 + self.filt_r]
				self.seq_range = self.wt_seq[self.start_pos: self.end_pos]
				break

			elif self.wt_seq[i:i+23] == rev_comp(self.t_seq):
				self.start_pos = i + 6 - self.end_range
				self.end_pos = i + 6 + self.end_range

				if self.start_pos < 0:
					s_num = (-1)*self.start_pos
					start_pos = 0

				if self.end_pos > len(self.wt_seq):
					e_num = self.end_pos - len(self.wt_seq)
					end_pos = len(self.wt_seq)

				self.s_seq = self.wt_seq[i + 6 - self.filt_r: i + 6 + self.filt_r]
				self.seq_range = self.wt_seq[self.start_pos: self.end_pos]
				break

		if self.start_pos == -1 or self.end_pos == -1:
			print("[Error] Can not find target sequence in reference sequence!!!")
			sys.exit()

		self.start_pos = self.start_pos
		self.end_pos = self.end_pos
		self.s_num = s_num
		self.e_num = e_num

		pri_for = self.seq_range[:15]
		pri_back = self.seq_range[-15:]

		self.for_table = [pri_for]
		self.back_table = [pri_back]

		for i in range(15):
			for n in 'ATGC':
				self.for_table.append(pri_for[:i] + n + pri_for[i+1:])
				self.back_table.append(pri_back[:i] + n + pri_back[i+1:])

	def read_file(self, obj):

		if obj == 'con':
			print('Read Control File...')
		elif obj == 'data':
			print('Read Data File...')

		self.range_len = len(self.seq_range)
		a=0
		b=0
		with open(self.input_file[obj]) as f:
			line_n = 0
			for line in f:
				line_n += 1
				if line_n%4 != 2:
					continue
				
				self.cnt_dict[obj]['all'] += 1
				
				if line.find('N') != -1: continue
				for_check = False
				back_check = False

				for pri in self.for_table:
					if line.find(pri) != -1:
						a+=1
						for_check = True
						start_seq = line.find(pri)
						break

				for back in self.back_table:
					if line.find(back) != -1:
						b+=1
						back_check = True
						end_seq = line.find(back) + 15
						break

				if for_check == True and back_check == True:
					self.cnt_dict[obj]['indi'] += 1
					line = line[start_seq:end_seq]
					if len(line) == self.range_len:
						self.cnt_dict[obj]['len'] += 1
						if line in self.line_dict[obj].keys():
							self.line_dict[obj][line] += 1
						else:
							self.line_dict[obj][line] = 1
					elif len(line) < self.range_len:
						self.cnt_dict[obj]['del'] += 1
					elif len(line) > self.range_len:
						self.cnt_dict[obj]['ins'] += 1
		print(a)
		print(b)

		value = {}
		for x in range(self.range_len + 1):
			value[x] = {'A':0, 'T':0, 'G':0, 'C':0}

		for seq, cnt in self.line_dict[obj].items():
			for order, nt in enumerate(seq):
				value[order][nt] += cnt

		tot = self.cnt_dict[obj]['len']

		rate_d = {}

		for x in range(self.range_len):
			ccv = 0
			for nt in 'ATGC':
				if nt != self.seq_range[x]:
					ccv += value[x][nt]
			rate_d[x] = ccv*100/tot

		self.rate_dict[obj] = rate_d
		self.value_dict[obj] = value

	def calculate_fold(self):

		self.rate_dict = {'data': [], 'con': []}

		for i in range(self.range_len):
			ccv = 0
			dcv = 0
			for x in 'ATGC':
				if x != self.seq_range[i]:
					ccv += self.value_dict['con'][i][x]
					dcv += self.value_dict['data'][i][x]
			self.rate_dict['con'].append(ccv*100/self.cnt_dict['con']['len'])
			self.rate_dict['data'].append(dcv*100/self.cnt_dict['data']['len'])

		self.fold_list = []

		for i in range(self.range_len):
			if self.rate_dict['con'][i] != 0:
				self.fold_list.append(self.rate_dict['data'][i]/self.rate_dict['con'][i])
			elif self.rate_dict['con'][i] == 0:
				self.fold_list.append(1.1)

	def curve_fitting(self):

		print("Fitting gaussian curve...")
		dict = {}

		for i in self.fold_list:
			i = int(i)
			if i in dict.keys():
				dict[i] += 1
			else:
				dict[i] = 1

		x = [-1, 0]
		y = [0,0]

		for a, b in dict.items():
			x.append(a+1)
			y.append(b)

		xn = [0]
		yn = [0]

		k = dict.keys()
		maxrange = max(x)
		for a in range(1, maxrange+1):
			if a in x:
				xn.append(a)
				yn.append(dict[a-1])
				continue
			xn.append(a)
			yn.append(0)

		r = 0
		while 0==0:
			try:
				popt, pcov = curve_fit(func, xn, yn, bounds=((-np.inf, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf)))
				break
			except:
				try:
					popt, pcov = curve_fit(func, xn, yn)
					break
				except:
					pass
			xn.append(maxrange + r)
			yn.append(0)
			r += 1

		a = popt[0]
		m = popt[1]
		c = popt[2]
		std = abs(popt[3])

		fold_n = m + 3.29 * std

		print('Mean : {0}'.format(m))
		print('std : {0}'.format(std))
		print('Fold : {0}'.format(round(fold_n,2)))

		self.fold = fold_n
		self.spot_list = []

		c = 0
		pos = -1
		for x in self.fold_list:
			pos += 1
			if x > fold_n:
				self.spot_list.append(pos)
				c += 1

		print(len(self.spot_list))

	def count_substitution(self):

		for seq, cnt in self.line_dict['data'].items():
			for s in self.spot_list:
				if self.seq_range[s] != seq[s]:
					self.cnt_dict['data']['sub'] += cnt
					break

	def fold_xlsx(self):

		print('make excel file ...')
		wb = xlsxwriter.Workbook(self.out_file+'_fold_result.xlsx')
		ws = wb.add_worksheet()

		ws.write('A1', 'control value')
		ws.write('A10', 'data value')
		ws.write('A19', 'data/control')
		ws.write('A22', 'control rate')
		ws.write('A30', 'data rate')

		row_list1 = []
		row_list2 = []
		for i in range(1, self.end_range + 1 - self.s_num):
			row_list1.append(i - self.end_range - 1 + self.s_num)

		for i in range(self.end_range + 1- self.s_num, self.end_range*2 + 1 -self.s_num - self.e_num):
			row_list2.append(i - self.end_range + self.s_num)

		for row in [0, 9, 21, 29, 17]:
			ws.write_row(row, 1, row_list1)
			ws.write_row(row, self.end_range + 1 - self.s_num, row_list2)

		col = 0
		for row in [1, 10, 22, 30]:
			for col in range(self.range_len):
				ws.write(row, col+1, self.seq_range[col])

		for col in [2, 11, 23, 31]:
			ws.write_column(0, col, ['A','T','G','C'])

		for i in range(self.range_len):
			row = 2
			ccv = 0
			dcv = 0
			for x in 'ATGC':
				ws.write(row, i+1, self.value_dict['con'][i][x])
				ws.write(row+9, i+1, self.value_dict['data'][i][x])
				if x != self.seq_range[i]:
					ccv += self.value_dict['con'][i][x]
					dcv += self.value_dict['data'][i][x]
				row += 1
			ws.write(row, i+1, round(ccv*100/self.cnt_dict['con']['len'],2))
			ws.write(row+9, i+1, round(dcv*100/self.cnt_dict['data']['len'],2))

		for i in range(self.range_len):
			row = 23
			for x in 'ATGC':
				ws.write(row, i+1, round(self.value_dict['con'][i][x]*100/self.cnt_dict['con']['len'],2))
				ws.write(row+9, i+1, round(self.value_dict['data'][i][x]*100/self.cnt_dict['data']['len'],2))
				row+=1

		for i in range(self.range_len):
			if self.rate_dict['con'][i] != 0:
				ws.write(18, i+1, round(self.rate_dict['data'][i]/self.rate_dict['con'][i], 2))
			elif self.rate_dict['con'][i] == 0:
				ws.write(18, i+1, 1.1)

		ws.write('I37', 'Fold')
		ws.write('I38', round(self.fold,2))
		ws.write_row('I39', list(self.cnt_dict['data'].keys()) + ['sub_per'])
		ws.write_row('I40', list(self.cnt_dict['data'].values()) + [round(self.cnt_dict['data']['sub']*100/self.cnt_dict['data']['indi'], 2)])

		print("Substitution rate: {0} %".format(round(self.cnt_dict['data']['sub']*100/self.cnt_dict['data']['indi'], 2)))

		chart = wb.add_chart({'type': 'scatter'})
		chart.add_series({'categories': ['Sheet1', 17, 1, 17, self.range_len], 'values': ['Sheet1', 18, 1, 18, self.range_len]})
		chart.set_title({'name': 'Data/Control'})
		ws.insert_chart('A37', chart)
		wb.close()


def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("wt_seq", type = str, help = "reference sequence(5' to 3')")
	parser.add_argument("target_seq", type = str, help = "target sequence with PAM (5' to 3')")
	parser.add_argument("data_file", type = str, help = "treated file for fastq-join format")
	parser.add_argument("con_file", type = str, help = "treated file for fastq-join format")
	parser.add_argument("out_file", type = str, help = "result file name")

	return parser.parse_args()

def main():

	tt = time.time()
	args = parse_args()

	n = NtSubstitution(args)

	n.seq_env()
	n.read_file('con')
	n.read_file('data')
	print(n.cnt_dict)
	n.calculate_fold()
	n.curve_fitting()
	n.count_substitution()
	n.fold_xlsx()

if __name__ == '__main__':
	main()
