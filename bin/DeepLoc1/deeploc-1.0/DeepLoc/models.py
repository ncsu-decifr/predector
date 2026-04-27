from theano import function as func
import theano.tensor as T
import lasagne
import pkg_resources
import numpy as np
from .utils import *

def net_tree():
	"""Compile a Convolutional BLSTM neural network for protein subcellular localization
	   
	Parameters:
		batch_size -- integer, minibatches size

	Outputs:
		deploy_fn -- compiled theano function for deploying
		l_out -- output of the network, can be used to save the model parameters		
	"""
	# Prepare Theano variables for inputs, masks and targets
	w_inits = lasagne.init.Orthogonal('relu')
	input_var = T.tensor3('inputs')
	mask_var = T.matrix('masks')

	n_feat = 23
	drop_per = 0
	drop_hid = 0
	n_filt = 20
	n_hid = 256
	n_class = 10
	# batch_size = 1

	# Input layer, holds the shape of the data
	l_in = lasagne.layers.InputLayer(shape=(None, None, n_feat), input_var=input_var)

	# Dropout positions of the protein sequence	
	l_indrop = DropoutSeqPosLayer(l_in, p=0)
	
	# Input layer with masks
	l_mask = lasagne.layers.InputLayer(shape=(None, None), input_var=mask_var)
	
	#### CNN ####
	# Size of convolutional layers
	f_size_a = 1
	f_size_b = 3
	f_size_c = 5
	f_size_d = 9
	f_size_e = 15
	f_size_f = 21
	
	# Shuffle shape to be properly read by the CNN layer
	l_shu = lasagne.layers.DimshuffleLayer(l_indrop, (0,2,1))
	
	# First convolutional layers
	l_conv_a = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1, W=w_inits, filter_size=f_size_a, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_b = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_b, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_c = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_c, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_d = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_d, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_e = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_e, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_f = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_f, nonlinearity=lasagne.nonlinearities.rectify)
	
	# Concatenate all CNN layers
	l_conc = lasagne.layers.ConcatLayer([l_conv_a, l_conv_b, l_conv_c, l_conv_d, l_conv_e, l_conv_f], axis=1)
	
	# Second CNN layer
	l_conv_final = lasagne.layers.Conv1DLayer(l_conc, num_filters=128, pad='same', stride=1, W=w_inits, filter_size=f_size_b, nonlinearity=lasagne.nonlinearities.rectify)
	
	# Reshuffle to initial shape
	l_indrop = lasagne.layers.DimshuffleLayer(l_conv_final, (0,2,1))
	
	# Dropout
	l_indrop = lasagne.layers.dropout(l_indrop, p=drop_hid)
	
	# LSTM forward and backward layers
	gate_parameters = lasagne.layers.recurrent.Gate(W_in=lasagne.init.Orthogonal(), W_hid=lasagne.init.Orthogonal(), 
		b=lasagne.init.Constant(0.))
	cell_parameters = lasagne.layers.recurrent.Gate(W_in=lasagne.init.Orthogonal(), W_hid=lasagne.init.Orthogonal(), W_cell=None,
		b=lasagne.init.Constant(0.),  nonlinearity=lasagne.nonlinearities.tanh)
	l_fwd = lasagne.layers.LSTMLayer(l_indrop, num_units=n_hid, name='LSTMFwd', mask_input=l_mask, ingate=gate_parameters, 
		forgetgate=gate_parameters, cell=cell_parameters, outgate=gate_parameters, nonlinearity=lasagne.nonlinearities.tanh, grad_clipping=2)
	l_bck = lasagne.layers.LSTMLayer(l_indrop, num_units=n_hid, name='LSTMBck', mask_input=l_mask, ingate=gate_parameters, 
		forgetgate=gate_parameters, cell=cell_parameters, outgate=gate_parameters, backwards=True, nonlinearity=lasagne.nonlinearities.tanh, 
		grad_clipping=2)
	
	# Concatenate both layers
	l_conc_lstm = lasagne.layers.ConcatLayer([l_fwd, l_bck], axis=2)
	
	# Attention mechanism
	l_dec = LSTMAttentionDecodeFeedbackLayer(l_conc_lstm, mask_input=l_mask, num_units=n_hid*2, aln_num_units=n_hid, n_decodesteps=10, name='LSTMAttention')
	
	# Slice last layer of the attention mechanism (context vector)
	l_last_hid = lasagne.layers.SliceLayer(l_dec, indices=-1, axis=1)

	# Final fully connected dense layer	
	l_dense = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_last_hid, p=drop_hid), name="Dense", num_units=n_hid*2, W=w_inits, nonlinearity=lasagne.nonlinearities.rectify)
	
	# Tree sigmoid layers
	l_sec = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Sec", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_ext = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Extracellular", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_int = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Internal", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_ergol = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="ER/Golgi", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_cmlys = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="CM/Lys", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_term = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Terminal", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_mitchl = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Mit/Chl", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_pernc = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Per/NC", nonlinearity=lasagne.nonlinearities.sigmoid)
	l_nuccyt = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Nuc/Cyt", nonlinearity=lasagne.nonlinearities.sigmoid)

	# Membrane output layer
	l_mem = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Sigmoid", W=w_inits, nonlinearity=lasagne.nonlinearities.sigmoid)	
	
	# Get output validation
	pt_sec, pt_ext, pt_int, pt_ergol, pt_cmlys, pt_term, pt_mitchl, pt_pernc, pt_nuccyt, test_membrane = lasagne.layers.get_output([l_sec, 
		l_ext, l_int, l_ergol, l_cmlys, l_term, l_mitchl, l_pernc, l_nuccyt, l_mem], inputs={l_in: input_var, l_mask: mask_var}, 
		deterministic=True)
	
	tree_test = T.concatenate([pt_sec.dimshuffle((0,'x')), pt_ext.dimshuffle((0,'x')), 
		pt_int.dimshuffle((0,'x')), pt_ergol.dimshuffle((0,'x')), pt_cmlys.dimshuffle((0,'x')),
		pt_term.dimshuffle((0,'x')), pt_mitchl.dimshuffle((0,'x')), pt_pernc.dimshuffle((0,'x')), 
		pt_nuccyt.dimshuffle((0,'x'))], axis=1)
	
	deploy_fn = func([input_var, mask_var], [tree_test, l_dec.alpha, test_membrane])
	return deploy_fn, [l_sec, l_ext, l_int, l_ergol, l_cmlys, l_term, l_mitchl, l_pernc, l_nuccyt, l_mem]

def net_softmax():
	"""Compile a Convolutional BLSTM neural network for protein subcellular localization
	   
	Parameters:
		batch_size -- integer, minibatches size

	Outputs:
		deploy_fn -- compiled theano function for deploying
		l_out -- output of the network, can be used to save the model parameters		
	"""
	# Prepare Theano variables for inputs, masks and targets
	w_inits = lasagne.init.Orthogonal('relu')
	input_var = T.tensor3('inputs')
	mask_var = T.matrix('masks')	

	n_feat = 23
	drop_per = 0
	drop_hid = 0
	n_filt = 20
	n_hid = 256
	n_class = 10
	# batch_size = 1

	# Input layer, holds the shape of the data
	l_in = lasagne.layers.InputLayer(shape=(None, None, n_feat), input_var=input_var)

	# Dropout positions of the protein sequence	
	l_indrop = DropoutSeqPosLayer(l_in, p=drop_per)
	
	# Input layer with masks
	l_mask = lasagne.layers.InputLayer(shape=(None, None), input_var=mask_var)
	
	#### CNN ####
	# Size of convolutional layers
	f_size_a = 1
	f_size_b = 3
	f_size_c = 5
	f_size_d = 9
	f_size_e = 15
	f_size_f = 21
	
	# Shuffle shape to be properly read by the CNN layer
	l_shu = lasagne.layers.DimshuffleLayer(l_indrop, (0,2,1))
	
	# First convolutional layers
	l_conv_a = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1, W=w_inits, filter_size=f_size_a, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_b = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_b, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_c = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_c, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_d = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_d, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_e = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_e, nonlinearity=lasagne.nonlinearities.rectify)
	l_conv_f = lasagne.layers.Conv1DLayer(l_shu, num_filters=n_filt, pad='same', stride=1,W=w_inits, filter_size=f_size_f, nonlinearity=lasagne.nonlinearities.rectify)
	
	# Concatenate all CNN layers
	l_conc = lasagne.layers.ConcatLayer([l_conv_a, l_conv_b, l_conv_c, l_conv_d, l_conv_e, l_conv_f], axis=1)
	
	# Second CNN layer
	l_conv_final = lasagne.layers.Conv1DLayer(l_conc, num_filters=128, pad='same', stride=1, W=w_inits, filter_size=f_size_b, nonlinearity=lasagne.nonlinearities.rectify)
	
	# Reshuffle to initial shape
	l_indrop = lasagne.layers.DimshuffleLayer(l_conv_final, (0,2,1))
	
	# Dropout
	l_indrop = lasagne.layers.dropout(l_indrop, p=drop_hid)
	
	# LSTM forward and backward layers
	gate_parameters = lasagne.layers.recurrent.Gate(W_in=lasagne.init.Orthogonal(), W_hid=lasagne.init.Orthogonal(), 
		b=lasagne.init.Constant(0.))
	cell_parameters = lasagne.layers.recurrent.Gate(W_in=lasagne.init.Orthogonal(), W_hid=lasagne.init.Orthogonal(), W_cell=None,
		b=lasagne.init.Constant(0.),  nonlinearity=lasagne.nonlinearities.tanh)

	l_fwd = lasagne.layers.LSTMLayer(l_indrop, num_units=n_hid, name='LSTMFwd', mask_input=l_mask, ingate=gate_parameters, 
		forgetgate=gate_parameters, cell=cell_parameters, outgate=gate_parameters, nonlinearity=lasagne.nonlinearities.tanh, grad_clipping=2)
	l_bck = lasagne.layers.LSTMLayer(l_indrop, num_units=n_hid, name='LSTMBck', mask_input=l_mask, ingate=gate_parameters, 
		forgetgate=gate_parameters, cell=cell_parameters, outgate=gate_parameters, backwards=True, nonlinearity=lasagne.nonlinearities.tanh,
		grad_clipping=2)
	
	# Concatenate both layers
	l_conc_lstm = lasagne.layers.ConcatLayer([l_fwd, l_bck], axis=2)
	
	# Attention mechanism
	l_dec = LSTMAttentionDecodeFeedbackLayer(l_conc_lstm, mask_input=l_mask, num_units=n_hid*2, aln_num_units=n_hid, n_decodesteps=10, name='LSTMAttention')
	
	# Slice last layer of the attention mechanism (context vector)
	l_last_hid = lasagne.layers.SliceLayer(l_dec, indices=-1, axis=1)

	# Final fully connected dense layer	
	l_dense = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_last_hid, p=drop_hid), name="Dense", num_units=n_hid*2, W=w_inits, nonlinearity=lasagne.nonlinearities.rectify)
	
	# Softmax output layer
	l_out = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=n_class, name="Softmax", W=w_inits, nonlinearity=lasagne.nonlinearities.softmax)

	# Membrane output layer
	l_mem = lasagne.layers.DenseLayer(lasagne.layers.dropout(l_dense, p=drop_hid), num_units=1, name="Sigmoid", W=w_inits, nonlinearity=lasagne.nonlinearities.sigmoid)	

	# Get output
	test_location, test_membrane = lasagne.layers.get_output([l_out, l_mem], 
		inputs={l_in: input_var, l_mask: mask_var}, deterministic=True)

	deploy_fn = func([input_var, mask_var], [test_location, l_dec.alpha, test_membrane])
	return deploy_fn, [l_out,l_mem]




def prediction(ids_array, protein_list, batch_size, attention=False):
	"""Function to perform a prediction on a input array
	   
	Parameters:
		ids_array -- list, protein ids
		protein_matrix -- np.array (float32), size: (N, maximum_length, 23), protein matrix in BLOSUM62 encoding
		mask_matrix -- np.array (float32), size: (N, maximum_length),
		batch_size -- integer, minibatches size
		attention -- bool, generate attention output

	Outputs:
		ids -- list, proteins ids (new order due to minibatch splitting)
		final_locations -- np.array (float32), predicted probability for the 10 possible localization
		final_membrane -- np.array (float32), predicted probability of a protein being membrane-bound or soluble
		final_alphas -- list, list with the attention values for each of the predicted proteins
	"""


	n_id = len(ids_array)
	n_proteins = len(protein_list)

	if not n_proteins == n_id:
		raise ValueError('The protein sequences, and ids have not the same number of examples')

	# Compile models
	deploy_fn1, network_out1 = net_tree()
	deploy_fn2, network_out2 = net_softmax()

	total_model = 16.0

	ids = np.array([])

	final_locations = np.array([]).reshape(0,10)
	final_membrane = np.array([])
	if attention:
		final_alphas = []


	# Iterate all the proteins	
	for batch in iterate_minibatches(ids_array, protein_list, batch_size):

		id, prot, mask = batch
		ids = np.concatenate((ids,id))
		n_prot = len(id)
		max_len = prot.shape[1]
		complete_tree = np.zeros((n_prot,9))
		complete_mem_test = np.zeros(n_prot)
		complete_soft = np.zeros((n_prot,10))
		if attention:
			complete_alpha = np.zeros((n_prot,max_len))

		# Itereate through the 16 models
		for e in [1,4]:
			for i in range(1,5):
				# Load parameters hierarchical tree
				params_net1 = pkg_resources.resource_filename('DeepLoc', 'parameters/params%s.%s.npz' % (i,e))
				with np.load(params_net1) as f1:
					param_values1 = [f1['arr_%d' % tt] for tt in range(len(f1.files))]

				lasagne.layers.set_all_param_values(network_out1, param_values1)

				# Load parameters softmax
				params_net2 = pkg_resources.resource_filename('DeepLoc', 'parameters/paramsx%s.%s.npz' % (i,e))
				with np.load(params_net2) as f2:
					param_values2 = [f2['arr_%d' % tt] for tt in range(len(f2.files))]

				lasagne.layers.set_all_param_values(network_out2, param_values2)

				# Predict the location probability, membrane-bound probability and attention values
				loc_tree, alpha_tree, membrane_tree = deploy_fn1(prot, mask)
				loc_soft, alpha_soft, membrane_soft = deploy_fn2(prot, mask)

				# Save all the values
				complete_tree += loc_tree
				complete_soft += loc_soft
				complete_mem_test += membrane_tree[0]
				complete_mem_test += membrane_soft[0]

				if attention:
					# Save attention values if needed
					complete_alpha += alpha_tree[0,-1,:]
					complete_alpha += alpha_soft[0,-1,:]

		# Average across the 16 models for the membrane-bound probability
		mem_pred = complete_mem_test / total_model
		mem_pred = np.around(mem_pred, decimals = 4)

		# Average for the tree and softmax model for the location probability
		pred_tree = complete_tree / (total_model / 2.0)
		pred_softmax = (tree_probs(pred_tree) + (complete_soft / (total_model / 2.0 ))) / 2.0

		# Reweghting based on class frequency (a priori information)
		reweight_softmax = Reweight(pred_softmax)

		final_locations = np.concatenate((final_locations, reweight_softmax), axis=0)
		final_membrane = np.concatenate((final_membrane, mem_pred), axis=0)
		
		if attention:
			alpha_pred = np.around(complete_alpha / total_model, decimals = 8)
			final_alphas.append(alpha_pred)

	if attention:
		return ids, final_locations, final_membrane, final_alphas
	else:
		return ids, final_locations, final_membrane
