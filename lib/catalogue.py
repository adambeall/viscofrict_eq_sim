# Code written by Dr Martijn van den Ende


import numpy as np
import matplotlib.pyplot as plt
from time import time

class generate_catalogue:

	V_c0 = 1e-2
	V_c1 = 1e-3

	def __init__(self, qdyn):
		self.q = qdyn
		pass

	def calc_magnitudes(self):
		self.get_events()
		self.calc_M0_catalogue()
		pass

	def calc_Mw(self, M0):
		return (2.0/3.0)*np.log10(M0) - 6.06

	def get_events(self):

		print("Catalogue: collecting events")

		ot = self.q.ot_vmax

		V0 = self.V_c0
		V1 = self.V_c1
		Vmax = ot["v"].values

		catalogue = []
		Nt = len(ot["t"])

		i0 = 0
		i1 = 1
		event = False

		"""
		Loop over time series
		If Vmax > V0: event starts
		If Vmax < V1: event stops
		Store (start, stop) indices of event
		"""

		for i in range(Nt):
			if Vmax[i] > V0 and event is False:
				i0 = i
				event = True
			if Vmax[i] < V1 and event is True:
				event = False
				i1 = i
				catalogue.append([i0, i1])

		self.catalogue = catalogue

	def calc_M0_catalogue(self):

		print("Catalogue: calculating moments")

		ot = self.q.ot_vmax
		ox = self.q.ox

		mask = np.isfinite(ox["x"])
		x = np.sort(ox["x"][mask].unique())
		t_vals = np.sort(ox["t"][mask].unique())
		t = ot["t"].values

		dx = x[1] - x[0]

		Nx = len(x)
		Nt = len(t_vals)-1

		v = ox["v"][mask][:Nx*Nt].values.reshape((Nt, Nx))
		vmax = ot["v"].values
		slip = ox["slip"][mask][:Nx*Nt].values.reshape((Nt, Nx))

		mu = self.q.set_dict["MU"]		

		catalogue = np.array(self.catalogue)
		M0 = np.zeros(len(catalogue))*np.nan
		duration = np.zeros(len(catalogue))*np.nan
		lengths = np.zeros(len(catalogue))*np.nan
		all_crack_inds = []

		"""
		Loop over all events in the catalogue
		For each event, find the fault elements that slipped with V(x) > threshold
		at some point during the event
		Compute average slip over the crack, and the length of the crack
		Assume circular crack to convert crack length to area
		"""

		for i, event in enumerate(catalogue):
			t0 = t[event[0]]
			t1 = t[event[1]]

			if t1 > t_vals[-1]: continue

			i0 = np.where(t_vals > t0)[0][0]
			i1 = np.where(t_vals > t1)[0][0]
			if i0 == i1: i1 += 1
			if i1 == len(slip): continue

			v_range = v[i0:i1]
			v_max = np.nanmax(v_range, axis=0)
			crack_inds = np.where(v_max > self.V_c0)[0]
			all_crack_inds.append(crack_inds)
			crack_slip = slip[i1][crack_inds] - slip[i0][crack_inds]
			avg_slip = np.mean(crack_slip)
			crack_length = len(crack_inds)*dx
			duration[i] = (t1-t0)/60.0
			M0[i] = mu*np.pi*avg_slip*(0.5*crack_length)**2
			lengths[i] = crack_length
		
		inds = np.isfinite(M0)
		self.finite_inds = inds
		self.M0 = np.array(M0[inds])
		self.lengths = np.array(lengths[inds])
		self.Mw = self.calc_Mw(self.M0)
		self.t = np.array([t[event[0]] for event in catalogue[inds]])
		self.crack_inds = all_crack_inds
		self.vmax = np.array([np.nanmax(vmax[event[0]:event[1]]) for event in catalogue[inds]])
		self.duration = np.array(duration[inds])
		del self.q
