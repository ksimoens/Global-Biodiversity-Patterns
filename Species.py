from parameters import*

class Species():

	def __init__(self,T,i):
		self.min_temp = T - NicheWidth
		self.max_temp = T + NicheWidth
		self.index = i
		self.alive = True
