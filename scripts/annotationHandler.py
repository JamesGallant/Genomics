import pandas as pd
import os
import argparse


def _fileReader(filename):
	file = os.path.join(os.getcwd(), filename)
	ext = os.path.splitext(file)[-1].lower()

	if ext == ".txt":
		return pd.read_table(file)
	elif ext == ".csv":
		return pd.read_csv(file)
	else:
		raise NameError("data must either be a .csv or a .txt file")


def _keepColumns(data):
	try:
		return pd.concat([data['Start'], data['Stop'], data['Locus tag'], data['Protein name']], axis=1)
	except Exception as e:
		return pd.concat([data['Start'], data['Stop'], data['Locus tag'], data['Protein Name']], axis=1)


def _editStrings(data, string_to_replace):
	try:
		data['Protein name'] = data['Protein name'].str.replace(string_to_replace, "_")
		return data
	except Exception as e:
		data['Protein name'] = data['Protein Name'].str.replace(string_to_replace, "_")
		data = data.drop('Protein Name', axis=1)
		return data


def _calculateGeneLength(data):
	geneLengthArr = []
	for start, stop in zip(data['Start'], data['Stop']):
		geneLength = stop-start
		geneLengthArr.append(geneLength)

	data.insert(4, 'Length', geneLengthArr)

	return data


def _addIntergenicRows(data):
	newdf = {'start': 		[],
			 'stop':  		[],
			 'gene':  		[],
			 'ORF':	  		[],
			 'Proteins':	[] }


	for idx, row in data.iterrows():
			newdf['start'].append(row['Start'])
			newdf['stop'].append(row['Stop'])
			newdf['gene'].append(row['Locus tag'])
			newdf['ORF'].append(row['ORF'])
			newdf['Proteins'].append(row['Protein name'])



			#add the intergenic space
			if idx < (len(data.index)-1):
				nextGeneStart = int(data.loc[idx + 1, 'Start'])
				intergenicStart = int(row['Stop'] + 1)
				intergenicStop = int(nextGeneStart - 1)

				if nextGeneStart > intergenicStart:
					newdf['start'].append(intergenicStart)
					newdf['stop'].append(intergenicStop)
					newdf['gene'].append('int')
					newdf['ORF'].append('INTERGENIC')
					newdf['Proteins'].append('n/a')

				else:
					continue
			else:
				break

	out = pd.DataFrame(newdf)		
	return out


def mainProcess(annotationFile, outDir):
	df = _fileReader(filename=annotationFile)
	df = _keepColumns(data=df)
	df = _editStrings(data=df, string_to_replace=" ")
	df = _editStrings(data=df, string_to_replace="/")
	df = _calculateGeneLength(data=df)
	df.insert(3, 'ORF', 'CDS')
	df = _addIntergenicRows(data=df)
	df.to_csv(outDir, sep="\t", header=False, index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Format NCBI annotation files for Pegasus")
	parser.add_argument('--file', required=True, help="Tabular annotation file of the reference")
	parser.add_argument('--out', required=True, help="output_directory")
	args = parser.parse_args()
	
	ANNOTATIONFILE = args.file
	OUTDIR = args.out

	
	mainProcess(annotationFile=ANNOTATIONFILE, outDir=OUTDIR)





