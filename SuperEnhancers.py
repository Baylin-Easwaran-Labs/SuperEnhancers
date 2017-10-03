import pandas as pd
import numpy as np
import re
import os
import glob
#from math import log
import multiprocessing
import warnings
#from sklearn.preprocessing import PolynomialFeatures
#from sklearn.metrics import r2_score
#from sklearn.metrics import mean_squared_error
#from sklearn.linear_model import LinearRegression
warnings.filterwarnings('ignore')

######################################################################
def MakeDirectory(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
######################################################################
def SaveDFtoCSV(directory,output_name, df):
	MakeDirectory(directory)
	path = os.path.join(directory, str(output_name))
	df.to_csv(path, sep=',', index=False)
######################################################################
def SaveDFtoPickle(directory, output_name, df):
	MakeDirectory(directory)
	path = os.path.join(directory, str(output_name))
	df.to_pickle(path)
######################################################################
def ListAllFiles(path,extension):
	file_array =[]
	for file in os.listdir(path):
		if file.endswith(extension):
			file_array.append(file)
	return file_array
######################################################################
def FinalDataFrameReconstruction(df_array):
    df_final = pd.DataFrame()
    for df1 in df_array:
        df_final = pd.concat([df_final, df1])
    print('Final DataFrame:\n{0}\n'.format(df_final.head()))
    return df_final
######################################################################
def MultiprocessingPreparationChrom(df):
	chrom_list = list(set(df['CHR']))
	df_array = []
	for chrom in chrom_list:
		df_t = df[df['CHR'] == chrom]
		if len(df_t) > 0:
			df_t1 = df_t.sort_values(by = ['Start'], ascending=True)
			df_array.append(df_t1)
	return chrom_list, df_array
######################################################################
def MultiprocsessingPerChrom(function,chrom_list, df_array, join_peak_distance):
	results = []
	results_async = []
	for chrom, df in zip(chrom_list, df_array):
		res = pool.apply_async(function, (chrom, df, join_peak_distance))
		results_async.append(res)

	results=[r.get() for r in results_async]

	return results
######################################################################
def MainWorker(chrom,df, join_peak_distance):
	#print(chrom)

	first = 0
	start_array = []
	end_array = []
	chrom_array = []
	peak_distance = []
	distance = 0
	for index, row in df.iterrows():
		start = row['Start']
		end=row['End']
		distance_temp = end - start
		if first == 0:
			start1 = start
			end1 = end
			distance = distance_temp
			first = 1
		elif first > 0:
			if (end1 + join_peak_distance) >= start:
				distance = distance+distance_temp
				end1 = end

			elif (end1 + join_peak_distance) < start:
				start_array.append(start1)
				end_array.append(end1)
				chrom_array.append(chrom)
				peak_distance.append(distance)

				start1 = start
				end1 = end
				distance = distance_temp

	df_final = pd.DataFrame({'CHR':chrom_array,'Start':start_array,'End':end_array,'Peaks_area':peak_distance})

	return df_final
######################################################################
def JoiningPeaks(folder, peak, join_peak_distance):
	
	path_peak = os.path.join(folder, peak)
	df_p = pd.read_csv(path_peak, sep ='\t', header = None)
	df_ps = df_p[[0,1,2]]
	df_ps.columns = ['CHR','Start','End']
	chrom_list, df_peaks_array = MultiprocessingPreparationChrom(df_ps)
	peaks_array = MultiprocsessingPerChrom(MainWorker,chrom_list, df_peaks_array, join_peak_distance)
	df_joined_peaks = FinalDataFrameReconstruction(peaks_array)
	print('Combining Peaks from file: {0}'.format(peak))
	print('Distance used for peak merging: {0}'.format(join_peak_distance))
	print('Initial Peak number: {0}'.format(len(df_p)))
	print('Merged Peaks number: {0}'.format(len(df_joined_peaks)))
	#name = peak +'_joined.csv'
	#SaveDFtoCSV('Results', name, df_joined_peaks)

	return df_joined_peaks
######################################################################
def CorrectingTiles(df, bedGraph_tile):
	#print('Correcting Tiles')
	array_start = []
	array_end = []
	array_CHR = []
	array_Signal = []
	
	for index, row in df.iterrows():

		chrom = row['CHR']
		start = row['Start']
		end = row['End']
		signal = row['Signal']
		while start < end:
			array_start.append(start)
			array_end.append(start+bedGraph_tile)
			array_CHR.append(chrom)
			array_Signal.append(signal)
			start = start + bedGraph_tile

	df_final = pd.DataFrame({'CHR':array_CHR,'Start':array_start,'End':array_end, 'Signal':array_Signal})
	return df_final
######################################################################
def TransformingBedGraphs(chrom, df, bedGraph_tile):
	print(chrom)
	
	df['Tile'] = df['Start']/bedGraph_tile
	##################################################
	max1 = (df['Tile'].max() - 1)
	if max1 > 0:
			pass
	else:
		max1 = 9970024

	array_tiles = list(range(0,int(max1)))
	array_chrom = [chrom]*int(max1)
	array_trick = [0]*int(max1)
	df_theo = pd.DataFrame({'CHR':array_chrom,'Tile':array_tiles,'Trick':array_trick})
	df_theo['Tile'] = df_theo['Tile'].astype(float)
	##################################################

	df_clean = df[df['Signal'] > 0]
	df_clean['Check'] = df_clean['End'] - df_clean['Start']
	df_clean_good = df_clean[df_clean['Check'] == bedGraph_tile]
	df_proc = df_clean[df_clean['Check'] != bedGraph_tile]
	
	if len(df_proc) > 0:
		df_cor = CorrectingTiles(df_proc, bedGraph_tile)
		df_clean_good = pd.concat([df_clean_good ,df_cor])

	try:
		df_final = df_clean_good[['CHR','Tile','Signal']]
		
		df_merge = pd.merge(df_theo,df_final, on = ['CHR','Tile'], how = 'outer')
		
		df_merge = df_merge[['CHR', 'Tile','Signal']]
	
		df_merge.fillna(0, inplace = True)
		
		df_merge.drop_duplicates(inplace = True)

		df_merge.sort_values(['Tile'], ascending=True)

		return df_merge
	except:
		return pd.DataFrame(columns=['CHR', 'Tile','Signal'])
######################################################################
def BedGraphManipulation(folder, bedgraph, bedGraph_tile):
	path = os.path.join(folder,bedgraph)
	print('\nReading Bedgraph File: {0}'.format(bedgraph))
	df = pd.read_csv(path, sep = '\t', names = ['CHR','Start','End','Signal'])

	#To be romeved
	df = df[df['CHR'] != 'chrM']
	#To be removed

	print('\nMultiprocessing Preparation')
	chrom_list, df_array = MultiprocessingPreparationChrom(df)
	print('\nStarting BedGraph Transformation')
	peaks_array = MultiprocsessingPerChrom(TransformingBedGraphs,chrom_list, df_array, bedGraph_tile)
	print('\nFormating DataFrame')
	df_final = FinalDataFrameReconstruction(peaks_array)

	return df_final	
######################################################################
def MultiprocessingPreparationChromTwoInput(df1,df2):
	chrom_list = list(set(df1['CHR']))
	df1_array = []
	df2_array = []
	for chrom in chrom_list:
		df_t1 = df1[df1['CHR'] == chrom]
		df_t2 = df2[df2['CHR'] == chrom]
		df1_array.append(df_t1)
		df2_array.append(df_t2)
	return chrom_list, df1_array,df2_array
######################################################################
def MultiprocsessingPerChromTwoInput(function,chrom_list, df1_array, df2_array):
	results = []
	results_async = []
	for chrom, df1, df2 in zip(chrom_list, df1_array, df2_array):
		res = pool.apply_async(function, (chrom, df1, df2))
		results_async.append(res)

	results=[r.get() for r in results_async]

	return results
######################################################################
def SignalMinusInput(df_s,df_in):
	chrom_list, df_s_array,df_i_array = MultiprocessingPreparationChromTwoInput(df_s,df_in)
	peaks_array = MultiprocsessingPerChromTwoInput(MainWorkerSignalMinusInput,chrom_list, df_s_array,df_i_array)
	df_final = FinalDataFrameReconstruction(peaks_array)

	return df_final
######################################################################
def MainWorkerSignalMinusInput(chrom, df_s,df_in):
	print('Merging Signal and Input on {0}'.format(chrom))
	if len(df_in) > 0:
		df_merge = pd.merge(df_s,df_in, on=['CHR','Tile'], how = 'outer')
		df_merge.fillna(0)
	
		print('Substracting Signal from Input on {0}'.format(chrom))
		df_merge['Diff'] = df_merge['Signal']-df_merge['Input']
		df_merge['Diff'][df_merge['Diff'] < 0] = 0
		df_merge = df_merge[['CHR','Tile','Diff']]
	else:
		df_merge = df_s[['CHR','Tile','Signal']]
		df_merge.columns['CHR','Tile','Diff']

	return df_merge
#####################################################################
def ReadsonPeaks(chrom, df_Peak, df_reads,bedGraph_tile):

	def ReadonPeaksMainWorker(start, end, peak_size, df_reads = df_reads, bedGraph_tile = bedGraph_tile):
		
		sum2 = 0
		sum3 = 0
		
		if len(df_reads) > 0:
			peak_size_with_empy_space = end-start

			s1 = start//bedGraph_tile
			e1 = end//bedGraph_tile

			s1_ = start/bedGraph_tile
			e1_ = end/bedGraph_tile
		
		
			df_read1 = df_reads[df_reads['Tile'] >= s1]
			df_read2 = df_read1[df_read1['Tile'] <= e1]
			df_read3 = df_read2[df_read2['Tile'] >= (s1+1)]

			if  s1 == s1_:
				df_read3 = df_read2

			if len(df_read3) == 0:
				print('THERE IS AN ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
	
			df_read3['Reads'] = df_read3['Diff']*bedGraph_tile
			sum1 = df_read3['Reads'].sum()

			start1 = df_read1['Diff'].iloc[0]
			end1 = df_read2['Diff'].iloc[-1]

			dif_s = ((s1+1) - s1_)*start1
			dif_end = (e1_- e1) * end1

			sum2 = (sum1+dif_s+dif_end)
			#sum3 = (sum1+dif_s+dif_end)/peak_size_with_empy_space
		
		return sum2

	df_Peak['Average_Reads'] = np.vectorize(ReadonPeaksMainWorker)(df_Peak['Start'],df_Peak['End'], df_Peak['Peaks_area'])
	
	return df_Peak
######################################################################
def MultiprocsessingPerChromTwoInputsStable(function,chrom_list, df1_array, df2_array, bedGraph_tile):
	results = []
	results_async = []
	for chrom, df1, df2 in zip(chrom_list, df1_array, df2_array):
		res = pool.apply_async(function, (chrom, df1, df2, bedGraph_tile))
		results_async.append(res)

	results=[r.get() for r in results_async]

	return results
######################################################################
def derivatives(x,x1, y,y1):
	run = x1-x
	rise = y1-y
	slope = rise/run
	return slope
######################################################################
def FindingXCoordinatesForSlope(x_array, y_array, slope1):
	i = 0
	slope_array = []
	while i < len(x_array):
		try:
			y = y_array[i]
			x = x_array[i]
			i = i + 1
			y1 = y_array[i]
			x1 = x_array[i]
			slope = derivatives(x,x1, y,y1)
			if slope >= slope1:
				slope_array.append(x[0])
		except:
			pass

	return slope_array
######################################################################
def CalculatingSuperEnhancers(df,slope,degree):

	Y = np.array(df['Reads_Norm'])
	Y = Y.reshape(-1, 1)
	X = np.array(df['Index'])
	X = X.reshape(-1, 1)
	pr = LinearRegression()
	quadratic = PolynomialFeatures(degree=degree)
	X_quad = quadratic.fit_transform(X)
	X_fit = np.arange(0,1,0.001)[:,np.newaxis]
	pr.fit(X_quad,Y)
	y_quad_fit=pr.predict(quadratic.fit_transform(X_fit))
	y_quad_pred = pr.predict(X_quad)
	print('Training MSE quadratic: % .3f' % (mean_squared_error(Y,y_quad_pred)))
	print('Training R^2 quadratic % .3f' %(r2_score(Y, y_quad_pred)))
	slope_array = FindingXCoordinatesForSlope(X_fit, y_quad_fit, slope)
	return slope_array
######################################################################
def FinilizeTheTables(df,column_name):
	df_final_peak = df[df[column_name] > 0]
	print('Finding the max reads')
	max_average_reads = df_final_peak[column_name].max()
	print('The max reads: {0}'.format(max_average_reads))
		
	print('Normilized Reads to 1')
	col1 = 'Reads_Norm_'+column_name
	df_final_peak[col1] = df_final_peak[column_name]/max_average_reads

	print('Sorting Enhancers according to Normalized Reads')
	df_final_peak.sort_values([col1], ascending=True, inplace = True)
	print('Reset Index')
	df_final_peak.reset_index(drop=True, inplace = True)
	print('Making New Index')
	df_final_peak['Index'] =df_final_peak.index
	print('Finding maximum number of Enhancers')
	max1 = df_final_peak['Index'].max()
	print('Normilize the number of Enhancers to 1')
	df_final_peak['Index'] = df_final_peak['Index']/max1
	print(df_final_peak.head())

	return df_final_peak
		
######################################################################
if __name__ == '__main__':
	join_peak_distance = 12500
	bedGraph_tile = 25
	#slope = 1
	#degree = 3

	pool = multiprocessing.Pool()

	folder_out = '/starter/starter-02/ikagiamp/data/Super_Enhancers/Data/Results_FINAL_NEW/'
	folder_bedGraph = '/starter/starter-02/ikagiamp/data/Super_Enhancers/Data/New_BedGraph_Galaxy/'
	folder_peakCall = '/starter/starter-02/ikagiamp/data/Super_Enhancers/Data/Peakcall/'
	
	
	array_peaks =\
	['Galaxy1743-[Mock_K27_Qnorm_bdgpeakcall_1500_300].txt',\
	'Galaxy1744-[455_K27_Qnorm_bdgpeakcall_1500_300].txt',\
	'Galaxy1745-[DAC_K27_Qnorm_bdgpeakcall_1500_300].txt',\
	'Galaxy1746-[DAC455_K27_Qnorm_bdgpeakcall_1500_300].txt']

	array_input_bed=\
	['Galaxy8-[Mock_INP_qNorm_Final].bedgraph',\
	'Galaxy5-[455_INP_qNorm_Final].bedgraph',\
	'Galaxy7-[DAC_INP_qNorm_Final].bedgraph',\
	'Galaxy6-[DAC455_INP_qNorm_Final].bedgraph']

	array_signal_bed=\
	['Galaxy12-[Mock_K27ac_Qnorm_Final].bedgraph',\
	'Galaxy9-[455_K27ac_Qnorm_Final].bedgraph',\
	'Galaxy11-[DAC_K27ac_Qnorm_Final].bedgraph',\
	'Galaxy10-[DAC455_K27ac_Qnorm_Final].bedgraph']
	
	
	#array_peaks =['Galaxy1744-[455_K27_Qnorm_bdgpeakcall_1500_300].txt']
	#array_input_bed=['Galaxy5-[455_INP_qNorm_Final].bedgraph']
	#array_signal_bed=['Galaxy9-[455_K27ac_Qnorm_Final].bedgraph']

	#array_peaks =['Galaxy1746-[DAC455_K27_Qnorm_bdgpeakcall_1500_300].txt']
	#array_input_bed=['Galaxy7-[DAC_INP_qNorm_Final].bedgraph']
	#array_signal_bed=['Galaxy10-[DAC455_K27ac_Qnorm_Final].bedgraph']

	for peak, bed_signal, bed_input in zip(array_peaks, array_signal_bed, array_input_bed):
		print(peak)
		print(bed_signal)
		print(bed_input)
		
		print('Joining peaks')
		peaks = JoiningPeaks(folder_peakCall, peak, join_peak_distance)
		chrom_set = list(set(peaks['CHR']))
		print('Number of Chromosomes after Peak Joining from the file {0} is: {1}'\
			.format(peak,len(chrom_set)))
		#print(chrom_set)
		#SaveDFtoCSV(folder_out, peak,peaks)
		
		
		print('Reading Signal')
		df_signal = BedGraphManipulation(folder_bedGraph, bed_signal, bedGraph_tile)
		print(df_signal.head())
		chrom_set = list(set(df_signal['CHR']))
		print('Number of Chromosomes after BedGraphManipulation Signal from the file {0} is: {1}'\
			.format(peak,len(chrom_set)))

		#name_10 = peak+'_pickle_signal'
		#SaveDFtoPickle(folder_out, name_10, df_signal)

		
		print('Reading Input')
		df_input = BedGraphManipulation(folder_bedGraph, bed_input, bedGraph_tile)
		print('Renaming Columns')
		df_input.columns = ['CHR','Tile','Input']
		print(df_input.head())
		chrom_set = list(set(df_input['CHR']))
		print('Number of Chromosomes after BedGraphManipulation Input from the file {0} is: {1}'\
			.format(peak,len(chrom_set)))

		#name_10 = peak+'_pickle_Input'
		#SaveDFtoPickle(folder_out, name_10, df_input)
		#break

		
		print('Signal Minus Input')
		df_c = SignalMinusInput(df_signal,df_input)
		df_c.fillna(0,inplace = True)
		chrom_set = list(set(df_c['CHR']))
		print('Number of Chromosomes after SignalMinusInput from the file {0} is: {1}'\
			.format(peak,len(chrom_set)))

		print('\nWorking on Peaks\n')
		print('Multiprocessing')
		chrom_list_peaks, df_array_peaks, df_array_c = MultiprocessingPreparationChromTwoInput(peaks,df_c)
		print('Reads on Peaks')
		df_array = MultiprocsessingPerChromTwoInputsStable(ReadsonPeaks,chrom_list_peaks, df_array_peaks, df_array_c, bedGraph_tile)
		print('Reconstructing DataFrame')
		df_final_peak = FinalDataFrameReconstruction(df_array)
		print('Working on Peaks Transformation')
		chrom_set = list(set(df_final_peak['CHR']))
		print('Number of Chromosomes after ReadsonPeaks from the file {0} is: {1}'\
			.format(peak,len(chrom_set)))

		df_peaks_Average_Reads = FinilizeTheTables(df_final_peak,'Average_Reads')
		name_enhancers1 = 'Enhancers_Average_Reads'+bed_signal+'.csv'
		chrom_set = list(set(df_peaks_Average_Reads['CHR']))
		print('Number of Chromosomes after FinilizeTheTables from the file {0} is: {1}'\
			.format(peak,len(chrom_set)))
		SaveDFtoCSV(folder_out, name_enhancers1,df_peaks_Average_Reads)

		#df_peaks_Average_Reads_With_Space = FinilizeTheTables(df_final_peak,'Average_Reads_With_Space')
		#name_enhancers2 = 'Enhancers_Average_Reads_With_Space'+bed_signal+'.csv'
		#SaveDFtoCSV(folder_out, name_enhancers2,df_peaks_Average_Reads_With_Space)
		
		
		#####################################################################
		'''
		print('Calculating Super Enhancers')
		slope_array = CalculatingSuperEnhancers(df_final_peak,slope,degree)

		if len(slope_array) > 0:
			limit = slope_array[0]
			print('Limit on the x-axis: {0}'.format(limit))
			df_super = df_final_peak[df_final_peak['Index'] >= limit]
			print('Number of Super Enhancers: {0}'.format(len(df_super)))
			df_enhancer = df_final_peak[df_final_peak['Index'] < limit]
			print('Number of Enhancers: {0}'.format(len(df_enhancer)))

			print('Saving')
			name_super = 'SuperEnhancers_'+bed_signal+'.csv'
			SaveDFtoCSV(folder_out, name_super,df_super)
			name_enhancers = 'Enhancers_'+bed_signal+'.csv'
			SaveDFtoCSV(folder_out, name_enhancers,df_enhancer)
		else:
			print('No super Enhancers')
			print('Number of Enhancers: {0}'.format(len(df_final_peak)))
			print('Saving')
			name_enhancers = 'Enhancers_'+bed_signal+'.csv'
			SaveDFtoCSV(folder_out, name_enhancers,df_final_peak)
		'''
		
