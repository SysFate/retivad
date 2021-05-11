'''
Created on 15 oct. 2020

@author: Francois STUDER

Variant detection:
@author: Stéfan

'''

import os
import tools
import regex
import argparse
from time import sleep,time
from threading import Thread,RLock
from Bio import SeqIO
from Bio.Seq import Seq
from sys import stderr as err
from multiprocessing import Pool
from datetime import datetime
samtools = 'samtools'
bedtools = 'bedtools'
bedops = 'bedops'
racon = 'racon'
minimap2 = 'minimap2'
medaka = 'medaka'

writeFile = RLock()
addInstance = RLock()
querrysGuibson = {}
sequenceGuibson = 'ACATTGAAGAACCTGTAGATAACTCGCTGT'
barcodes5 = {}
barcodes3 = {}
statsLength = {'length_between_Guibson' : {}, 'length_original_seq': {}, 'length_tooShort_seq': {}, 'length_tooLong_seq': {}, 'length_match_align': {}, 'length_good_insert': {}}
statsLongueurNBGuib = [[],[]]
stats = {}
statsShort = {'Seq before cut' : 0, 'Seq without guibson' : 0, 'Seq too short' : 0, 'Seq too long' : 0, 'Seq without barcode cDNA' : 0, 'Multiple barcode cDNA' : 0, 'Seq without barcode PCR' : 0, 'Multiple barcode PCR' : 0, 'Seq with issues in quality control' : 0, 'Total good sequences' : 0,'Total Seq' : 0}
statsBarcodesGood = {}
extractedBarcode = set()
writeStats = RLock()
alignment = RLock()
statsAlignementLock = RLock()
statsAlignement = {}

def find_cut_sequence(record,querrys,lengthMini):
    listExtract = []
    listExtractMatch = []
    find_cut_sequence_rename(record,querrys,listExtract,listExtractMatch,lengthMini)
    return listExtract,listExtractMatch

def find_cut_sequence_rename(record,querrys,listExtract,listExtractMatch,lengthMini,realBeginPosition=0,beginSearch=[]):
    if len(record) >= lengthMini:
        find,beginSearch = find_best_match_entry(record,querrys,beginSearch=beginSearch)
        if find != None:
            begin = 0
            for coordinate in find.spans():
                find_cut_sequence_rename(record[begin:coordinate[0]],querrys,listExtract,listExtractMatch,lengthMini,realBeginPosition=realBeginPosition+begin,beginSearch=beginSearch)
                listExtractMatch.append(record[coordinate[0]:coordinate[1]])
                listExtractMatch[-1].id = listExtractMatch[-1].id+'_'+str(realBeginPosition+coordinate[0])+'_'+str(realBeginPosition+coordinate[1])
                begin = coordinate[1]
            find_cut_sequence_rename(record[begin:],querrys,listExtract,listExtractMatch,lengthMini,realBeginPosition=realBeginPosition+begin,beginSearch=beginSearch)
        else:
            record = record[:]
            record.id = record.id+'_'+str(realBeginPosition)+'_'+str(realBeginPosition+len(record))
            listExtract.append(record)
    else:
        record = record[:]
        record.id = record.id+'_'+str(realBeginPosition)+'_'+str(realBeginPosition+len(record))
        listExtract.append(record)
        
def find_best_match_entry(record,querrys,beginSearch=[]):
    find = regex.search(querrys[-1][0], str(record.seq),concurrent=True)
    if find != None and len(find.spans()) != 0:
        if len(querrys) > 1 and (len(beginSearch) == 0 or beginSearch[0] != (len(querrys)-1)):
            newFind,beginSearchSub = find_best_match(record,querrys[:-1],begin=beginSearch)
            if newFind != None:
                return newFind,beginSearchSub
            elif len(querrys[-1][1]) > 0:
                newFind,beginSearchSub = find_best_match(record,querrys[-1][1],begin=beginSearch[1:])
                if newFind != None:
                    newBegin = [len(querrys)-1]
                    newBegin.extend(beginSearchSub)
                    return newFind,newBegin
                else:
                    return find,[len(querrys)-1,-1]
            else:
                return find,[len(querrys)-1]
        elif len(querrys[-1][1]) > 0:
            newFind,beginSearchSub = find_best_match(record,querrys[-1][1],begin=beginSearch[1:])
            if newFind != None:
                newBegin = [len(querrys)-1]
                newBegin.extend(beginSearchSub)
                return newFind,newBegin
            else:
                return find,[len(querrys)-1,-1]
        else:
            return find,[len(querrys)-1]
    else:
        return None,[-1]

def find_best_match(record,querrys,begin=None):
    # On cherche pour la meilleur valeur de la list les valeurs sont toute decroissante
    if begin != None and len(begin) > 0:
        if begin[0] == -1:
            return None,begin
        i = begin[0]
    else:
        i = 0
    find = regex.search(querrys[i][0], str(record.seq),concurrent=True)
    while (find == None or len(find.spans()) == 0) and len(querrys) > i+1:
        i = i+1
        find = regex.search(querrys[i][0], str(record.seq),concurrent=True)
    if find != None and len(find.spans()) != 0:
        if len(begin) > 0 and i != begin[0]:
            begin = []
        newbegin = [i]
        if len(querrys[i][1]) > 0:
            newFind,newBeginUnd = find_best_match(record,querrys[i][1],begin=begin[1:])
            newbegin.extend(newBeginUnd)
            if newFind != None:
                return newFind,newbegin
            else:
                return find,newbegin
        else:
            return find,newbegin
    else:
        return None,[-1]

def search_sequences(record,listExtract,querry,realpos=0,secondSearch=False):
    find = regex.search(querry, str(record.seq),concurrent=True)
    if find != None and len(find.spans()) != 0:
        listExtract.append([find.spans()[0][0]+realpos,find.spans()[0][1]+realpos])
        if secondSearch:
            return True
        else:
            #if search_sequences(record[:find.spans()[0][0]],listExtract,querry,secondSearch=True):
            #    return True
            #else:
            return search_sequences(record[find.spans()[0][1]:],listExtract,querry,realpos=find.spans()[0][1]+realpos,secondSearch=True)
    else:
        return False
    
def search_barcode(record,barcodes):
    listBarcodeDetected = []
    barcordesExtracts = {}
    numberBarcodeDetected = 0
    for barcodeName,barcodeQuerryStrand in barcodes.items():
        barcordesExtracts[barcodeName] = {}
        for strand,querry in barcodeQuerryStrand.items():
            barcordesExtracts[barcodeName][strand] = []
            search_sequences(record,barcordesExtracts[barcodeName][strand],querry)
            if len(barcordesExtracts[barcodeName][strand]) > 0:
                numberBarcodeDetected = numberBarcodeDetected + len(barcordesExtracts[barcodeName][strand])
                listBarcodeDetected.append(barcodeName+strand)
                if numberBarcodeDetected > 1:
                    return listBarcodeDetected,barcordesExtracts,numberBarcodeDetected
            else:
                del barcordesExtracts[barcodeName][strand]
        if len(barcordesExtracts[barcodeName]) == 0:
            del barcordesExtracts[barcodeName]
    
    return listBarcodeDetected,barcordesExtracts,numberBarcodeDetected

def update_record_id(record,dicoExtracts):
    record.id = record.id+"_("
    for barcodeName,strandData in dicoExtracts.items():
        for strand,positions in strandData.items():
            for pos in positions:
                record.id = record.id+"%"+barcodeName+strand+"\\"+str(pos[0])+","+str(pos[1])
            
    record.id = record.id+")"

def cut_guibson(record,querrysGuibson,lengthMiniGuibson):
    resultsAfterCut = []
    statsSeq = {'length':len(record),'posGuib':[]}
    mapGuibSeq = [False]*len(record)
    listExtract,resultsSeqGuib = find_cut_sequence(record,querrysGuibson['+'],lengthMiniGuibson)
    for guibSeq in resultsSeqGuib:
        datasGuib = guibSeq.id.split("_")
        mapGuibSeq[int(datasGuib[1]):int(datasGuib[2])] = [True]*(int(datasGuib[2])-int(datasGuib[1]))
        statsSeq['posGuib'].append([int(datasGuib[1]),int(datasGuib[2]),float(datasGuib[1])/len(record)*100,float(datasGuib[2])/len(record)*100])
    for extract in listExtract:
        listExtractCut,listExtractMatch = find_cut_sequence(extract,querrysGuibson['-'],lengthMiniGuibson)
        resultsAfterCut.extend(listExtractCut)
        resultsSeqGuib.extend(listExtractMatch)
        for guibSeq in listExtractMatch:
            datasGuib = guibSeq.id.split("_")
            mapGuibSeq[int(datasGuib[1])+int(datasGuib[3]):int(datasGuib[1])+int(datasGuib[4])] = [True]*(int(datasGuib[4])-int(datasGuib[3]))
            statsSeq['posGuib'].append([int(datasGuib[3])+int(datasGuib[1]),int(datasGuib[4])+int(datasGuib[1]),(float(datasGuib[3])+float(datasGuib[1]))/len(record)*100,(float(datasGuib[4])+float(datasGuib[1]))/len(record)*100])
    # Ici on veut faire des stats sur la distance entre les guibsons.
    lengthBetGuib = {}
    if len(resultsSeqGuib) > 1:
        firstGuib = True
        counting = False
        for case in mapGuibSeq:
            if case:
                if firstGuib:
                    firstGuib = False
                    counting = False
                elif counting:
                    if i in lengthBetGuib:
                        lengthBetGuib[i] = lengthBetGuib[i] + 1
                    else:
                        lengthBetGuib[i] = 1
                    counting = False
            elif counting:
                i = i+1
            else:
                i = 1
                counting = True

    statsSeq['nbGuib'] = len(resultsSeqGuib)
                
    return resultsAfterCut,statsSeq,lengthBetGuib,record.id

def test_sequence(record,barcodes5,barcodes3,miniLength,maxLength,distanceMini,distanceMax):
    if len(record) < miniLength:
        lengthFile = -1
        nameBarcode = ''
        numberBarcode5Detected = 0
        numberBarcode3Detected = 0
        typeMatch = 0
    elif len(record) > maxLength:
        lengthFile = 1
        nameBarcode = ''
        numberBarcode5Detected = 0
        numberBarcode3Detected = 0
        typeMatch = 0
    else:
        lengthFile = 0
        listBarcode5Detected,barcordes5Extracts,numberBarcode5Detected = search_barcode(record,barcodes5)
        
        if numberBarcode5Detected == 1:
            listBarcode3Detected,barcordes3Extracts,numberBarcode3Detected = search_barcode(record,barcodes3)
            
            if numberBarcode3Detected == 1:
                nameBarcode = listBarcode5Detected[0][:-1]+"_"+listBarcode3Detected[0][:-1]
                update_record_id(record,barcordes5Extracts)
                update_record_id(record,barcordes3Extracts)
                if listBarcode5Detected[0][-1] == listBarcode3Detected[0][-1]:
                    positionsBC5 = barcordes5Extracts[listBarcode5Detected[0][:-1]][listBarcode5Detected[0][-1]][0]
                    positionsBC3 = barcordes3Extracts[listBarcode3Detected[0][:-1]][listBarcode3Detected[0][-1]][0]
                    if listBarcode5Detected[0][-1] == "+":
                        distance = positionsBC3[0]-positionsBC5[1]
                    else:
                        distance = positionsBC5[0]-positionsBC3[1]
                    if distance < 0:
                        typeMatch = 5
                    elif distance < distanceMini:
                        typeMatch = 2
                    elif distance > distanceMax:
                        typeMatch = 3
                    else:
                        typeMatch = 1
                else:
                    typeMatch = 4
            else:
                nameBarcode = ''
                typeMatch = 0
        else:
            nameBarcode = ''
            numberBarcode3Detected = -1
            typeMatch = 0
    
    return record,nameBarcode,numberBarcode5Detected,numberBarcode3Detected,lengthFile,typeMatch

class annalyseRecord(Thread):
    '''
    classdocs
    '''


    def __init__(self, listRecord,outFolder):
        '''
        Constructor
        '''
        Thread.__init__(self)
        self.listRecord = listRecord
        self.outComplex = outFolder
        
    def run(self):
        global calculPool
        global lengthMiniGuibson
        global querrysGuibson
        global barcodes5
        global barcodes3
        global miniLength
        global maxLength
        global distanceMini
        global distanceMax
        global stats
        global statsShort
        global writeFile
        global writeStats
        global extractedBarcode
        global statsLength
        listJob = []
        localStatsLength = {'length_between_Guibson' : {}, 'length_original_seq': {}, 'length_tooShort_seq': {}, 'length_tooLong_seq': {}, 'length_match_align': {}, 'length_good_insert': {}}
        for record in self.listRecord:
            job = calculPool.apply_async(cut_guibson, (record,querrysGuibson,lengthMiniGuibson))
            listJob.append(job)
        
        localStatsBySeq = {}
        numberSeqAnnalyse = len(self.listRecord)
        self.listRecord = None
        results = {'withoutGuib' : [], 'multi' : {5 : [], 3 : []}, 'noGoodLength' : {-1 : [], 1 : []}, 'noBarcode' : {5 : [], 3 : []}, 'Alone 3\'' : {-1 : [], 0 : [], 1 : []}}
        listSecondJob = []
        for job in listJob:
            result = job.get()
            if result[1]['nbGuib'] == 0:
                results['withoutGuib'].append(result[0][0].format('fastq'))
            else:
                for record in result[0]:
                    job = calculPool.apply_async(test_sequence, (record,barcodes5,barcodes3,miniLength,maxLength,distanceMini,distanceMax))
                    listSecondJob.append(job)
                for length,number in result[2].items():
                    if length in localStatsLength['length_between_Guibson']:
                        localStatsLength['length_between_Guibson'][length] += number
                    else:
                        localStatsLength['length_between_Guibson'][length] = number
            localStatsBySeq[result[3]] = result[1]
            if result[1]['length'] in localStatsLength['length_original_seq']:
                localStatsLength['length_original_seq'][result[1]['length']] += 1
            else:
                localStatsLength['length_original_seq'][result[1]['length']] = 1
        
        listJob = None
        resultsBarcode = {}
        for nameBarcode5,querry in barcodes5.items():
            for nameBarcode3,querry in barcodes3.items():
                resultsBarcode[nameBarcode5+"_"+nameBarcode3] = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 'uncutRecord': [], 'uncutRecordBugStrand': [], 'uncutRecordBugCoordonne': [], 'goodInsert': [], 'trashSeq': []}
        for job in listSecondJob:
            result = job.get()
            if result[4] == 0:
                if result[2] == 0:
                    results['noBarcode'][5].append(result[0].format('fastq'))
                elif result[2] == 1:
                    if result[3] == 0:
                        results['noBarcode'][3].append(result[0].format('fastq'))
                    elif result[3] == 1:
                        resultsBarcode[result[1]][result[5]].append(result[0].format('fastq'))
                        if result[5] == 1:
                            record = result[0]
                            idCut = record.id.split('(')
                            if len(idCut) == 3:
                                posBC5 = idCut[1].split("%")
                                if len(posBC5) == 2 and posBC5[0] == '':
                                    BC5Coup = posBC5[1].split("\\")
                                    posBC5 = BC5Coup[1].split(',')
                                else:
                                    posBC5 = []
                                posBC3 = idCut[2].split("%")
                                if len(posBC3) == 2 and posBC3[0] == '':
                                    BC3Coup = posBC3[1].split("\\")
                                    posBC3 = BC3Coup[1].split(',')
                                else:
                                    posBC3 = []
                                if len(posBC3) == 0 or len(posBC5) == 0:
                                    resultsBarcode[result[1]]['uncutRecordBugCoordonne'].append(record.format('fastq'))
                                elif BC5Coup[0][-1] == BC3Coup[0][-1]:
                                    if BC5Coup[0][-1] == '+':
                                        if posBC5[0] != '0':
                                            recordCut = record[:int(posBC5[0])]
                                            recordCut.id = recordCut.id+"_0,"+posBC5[0]
                                            resultsBarcode[result[1]]['trashSeq'].append(recordCut.format('fastq'))
                                        goodInsert = record[int(posBC5[1].split(')')[0]):int(posBC3[0])].format('fastq')
                                        resultsBarcode[result[1]]['goodInsert'].append(goodInsert)
                                        posEnd = int(posBC3[1].split(')')[0])
                                        if posEnd != len(record):
                                            recordCut = record[posEnd:]
                                            recordCut.id = recordCut.id+"_"+str(posEnd)+","+str(len(record))
                                            resultsBarcode[result[1]]['trashSeq'].append(recordCut.format('fastq'))
                                    elif BC3Coup[0][-1] == '-':
                                        if posBC3[0] != '0':
                                            recordCut = record[:int(posBC3[0])]
                                            recordCut.id = recordCut.id+"_0,"+posBC3[0]
                                            resultsBarcode[result[1]]['trashSeq'].append(recordCut.format('fastq'))
                                        goodInsert = record[int(posBC3[1].split(')')[0]):int(posBC5[0])].format('fastq')
                                        resultsBarcode[result[1]]['goodInsert'].append(goodInsert)
                                        posEnd = int(posBC5[1].split(')')[0])
                                        if posEnd != len(record):
                                            recordCut = record[posEnd:]
                                            recordCut.id = recordCut.id+"_"+str(posEnd)+","+str(len(record))
                                            resultsBarcode[result[1]]['trashSeq'].append(recordCut.format('fastq'))
                                    else:
                                        resultsBarcode[result[1]]['uncutRecordBugStrand'].append(record.format('fastq'))
                                        goodInsert = None
                                    if goodInsert != None:
                                        if len(goodInsert) in localStatsLength['length_good_insert']:
                                            localStatsLength['length_good_insert'][len(goodInsert)] += 1
                                        else:
                                            localStatsLength['length_good_insert'][len(goodInsert)] = 1
                                else:
                                    resultsBarcode[result[1]]['uncutRecordBugStrand'].append(record.format('fastq'))
                            else:
                                resultsBarcode[result[1]]['uncutRecord'].append(record.format('fastq'))
                    else:
                        results['multi'][3].append(result[0].format('fastq'))
                else:
                        results['multi'][5].append(result[0].format('fastq'))
            else:
                results['noGoodLength'][result[4]].append(result[0].format('fastq'))
                if result[4] == -1:
                    if len(result[0]) in localStatsLength['length_tooShort_seq']:
                        localStatsLength['length_tooShort_seq'][len(result[0])] += 1
                    else:
                        localStatsLength['length_tooShort_seq'][len(result[0])] = 1
                elif result[4] == 1:
                    if len(result[0]) in localStatsLength['length_tooLong_seq']:
                        localStatsLength['length_tooLong_seq'][len(result[0])] += 1
                    else:
                        localStatsLength['length_tooLong_seq'][len(result[0])] = 1
        
        seqStatGeneralToWrite = []
        seqStatPosGuibToWrite = []
        localStatsLongueurNBGuib = [[],[]]
        if args.dev:
            for seqName,datas in localStatsBySeq.items():
                seqStatGeneralToWrite.append('{}\t{}\t{}'.format(seqName,str(datas['length']),str(datas['nbGuib'])))
                localStatsLongueurNBGuib[0].append(datas['length'])
                localStatsLongueurNBGuib[1].append(datas['nbGuib'])
                for positions in datas['posGuib']:
                    seqStatPosGuibToWrite.append('{}\t{}\t{}\t{}\t{}'.format(seqName,str(positions[0]),str(positions[1]),str(positions[2]),str(positions[3])))
        else:
            for seqName,datas in localStatsBySeq.items():
                seqStatGeneralToWrite.append('{}\t{}\t{}'.format(seqName,str(datas['length']),str(datas['nbGuib'])))
                localStatsLongueurNBGuib[0].append(datas['length'])
                localStatsLongueurNBGuib[1].append(datas['nbGuib'])
            
        with writeFile:
            if args.dev:
                for nameBarcode,datas in resultsBarcode.items():
                    if len(datas[0]) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_unknowGood.fastq'),datas[0])
                    if len(datas[1]) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'.fastq'),datas[1])
                    if len(datas[2]) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_tooClose.fastq'),datas[2])
                    if len(datas[3]) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_tooFar.fastq'),datas[3])
                    if len(datas[4]) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_strandProb.fastq'),datas[4])
                    if len(datas[5]) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_positionProb.fastq'),datas[5])
                    if len(datas['uncutRecord']) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_uncut.fastq'),datas['uncutRecord'])
                    if len(datas['uncutRecordBugStrand']) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_uncutStrandProb.fastq'),datas['uncutRecordBugStrand'])
                    if len(datas['uncutRecordBugCoordonne']) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_uncutCoordonne.fastq'),datas['uncutRecordBugCoordonne'])
                    if len(datas['trashSeq']) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_trash.fastq'),datas['trashSeq'])
                    if len(datas['goodInsert']) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_cut.fastq'),datas['goodInsert'])
                        extractedBarcode.add(nameBarcode)

                if len(results['multi'][5]) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'multiBarcode_5.fastq'),results['multi'][5])
                if len(results['multi'][3]) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'multiBarcode_3.fastq'),results['multi'][3])
                
                if len(results['noGoodLength'][-1]) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'littleSeq.fastq'),results['noGoodLength'][-1])
                if len(results['noGoodLength'][1]) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'longSeq.fastq'),results['noGoodLength'][1])
                
                if len(results['noBarcode'][5]) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'withoutBarcode_5.fastq'),results['noBarcode'][5])
                if len(results['noBarcode'][3]) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'withoutBarcode_3.fastq'),results['noBarcode'][3])
                    
                if len(results['withoutGuib']) != 0:
                    tools.append_File_text(os.path.join(self.outComplex,'withoutGuibson.fastq'),results['withoutGuib'])
                tools.append_File_text(os.path.join(self.outComplex,'statOriginalSeqPosGuib.tsv'),seqStatPosGuibToWrite)
            else:
                for nameBarcode,datas in resultsBarcode.items():
                    if len(datas['goodInsert']) != 0:
                        tools.append_File_text(os.path.join(self.outComplex,nameBarcode+'_cut.fastq'),datas['goodInsert'])
                        extractedBarcode.add(nameBarcode)
            tools.append_File_text(os.path.join(self.outComplex,'statOriginalSeq.tsv'),seqStatGeneralToWrite)
        
        with writeStats:
            for nameBarcode,datas in resultsBarcode.items():
                stats['Number good unknow '+nameBarcode] += len(datas[0])
                stats['Number good '+nameBarcode] += len(datas[1])
                stats['Number good too close '+nameBarcode] += len(datas[2])
                stats['Number good too far '+nameBarcode] += len(datas[3])
                stats['Number good bad strand '+nameBarcode] += len(datas[4])
                stats['Number good bad position '+nameBarcode] += len(datas[5])
                
                statsShort['Total good sequences'] += len(datas[1])
                statsBarcodesGood['Number good '+nameBarcode] += len(datas[1])
                statsShort['Seq with issues in quality control'] = statsShort['Seq with issues in quality control'] + len(datas[0]) + len(datas[2]) + len(datas[3]) + len(datas[4]) + len(datas[5])
                statsShort['Total Seq'] = statsShort['Total Seq'] + len(datas[0]) + len(datas[1]) + len(datas[2]) + len(datas[3]) + len(datas[4]) + len(datas[5])
            
            statsShort['Seq too short'] += len(results['noGoodLength'][-1])
            statsShort['Seq too long'] += len(results['noGoodLength'][1])

            statsShort['Multiple barcode cDNA'] += len(results['multi'][5])
            statsShort['Multiple barcode PCR'] += len(results['multi'][3])
            statsShort['Seq without barcode cDNA'] += len(results['noBarcode'][5])
            statsShort['Seq without barcode PCR'] += len(results['noBarcode'][3])
            statsShort['Total Seq'] = statsShort['Total Seq'] + len(results['noBarcode'][5]) + len(results['noBarcode'][3]) + len(results['noGoodLength'][1]) + len(results['noGoodLength'][-1]) + len(results['multi'][5]) + len(results['multi'][3]) + len(results['withoutGuib'])
                
            statsShort['Seq before cut'] = statsShort['Seq before cut'] + numberSeqAnnalyse
            statsShort['Seq without guibson'] = statsShort['Seq without guibson'] + len(results['withoutGuib'])
            
            for nameStats,datas in localStatsLength.items():
                for length,number in datas.items():
                    if length in statsLength[nameStats]:
                        statsLength[nameStats][length] += number
                    else:
                        statsLength[nameStats][length] = number
            
            statsLongueurNBGuib[0].extend(localStatsLongueurNBGuib[0])
            statsLongueurNBGuib[1].extend(localStatsLongueurNBGuib[1])
            
        resultsBarcode = None
        listSecondJob = None
        results = None
        record = None
        recordCut = None
        localStatsLength = None
        goodInsert = None
        localStatsBySeq = None
        seqStatGeneralToWrite = None
        seqStatPosGuibToWrite = None
        localStatsLongueurNBGuib = None
        self.outComplex = None

def annalyse_file(fastqToAnnalyse,outFolder):
    global listInstance
    global core
    listLocalInstance = []
    listRecord = []
    for record in SeqIO.parse(fastqToAnnalyse, "fastq"):
        listRecord.append(record)
        if len(listRecord) == numberRecordByAnnalyse:
            listLocalInstance.append(annalyseRecord(listRecord,outFolder))
            with addInstance:
                tools.auto_add_thread(listInstance,0.5,listLocalInstance[-1])
            listRecord = []
            sleep(0.1)
    if len(listRecord) != 0:
        listLocalInstance.append(annalyseRecord(listRecord,outFolder))
        with addInstance:
            tools.auto_add_thread(listInstance,1,listLocalInstance[-1])
    tools.wait_end_thread_list(listLocalInstance)

def generate_aligner_command():
    cmdGenerateBed = "cat fileAlign"
    cmdGenerateBed = " ".join([cmdGenerateBed, "|", samtools, "view", "-@", str(core), "-Sb1", "-", "|", samtools, "sort", "-T", os.path.join("outFolder","barcodeName.samsort"), "-@", str(core), "-"])
    cmdGenerateBed = " ".join([cmdGenerateBed, "|", bedtools, "bamtobed", "-i", "-"])
    cmdGenerateBed = " ".join([cmdGenerateBed, ">", os.path.join("outFolder","barcodeName.bed")])
    if args.dev:
        tools.make_dirs(os.path.join(devFolder,"bedgraph"))
        cmdGenerateBedgraph = " ".join([bedtools,"genomecov -bga -split -i",os.path.join("outFolder","barcodeName.bed"),"-g",args.genome+'.chromsize','>',os.path.join(devFolder,"bedgraph","barcodeName.bedgraph")])
    else:
        cmdGenerateBedgraph = " ".join([bedtools,"genomecov -bga -split -i",os.path.join("outFolder","barcodeName.bed"),"-g",args.genome+'.chromsize','>',os.path.join("outFolder","barcodeName.bedgraph")])

    if args.aligner == 'bwa':
        # mem -k 10 -W 10 -r 5 -A 1 -B 1 -O 1 -E 1 -L 0 -y 14 -c 1000 -m 100 -t 80
        if args.alignerSpecialsParam == None:
            specialsParam = "-x ont2d"
        else:
            specialsParam = args.alignerSpecialsParam
        if args.alignerPath == None:
            alignerPath = args.aligner
        else:
            alignerPath = args.alignerPath
        return " ".join([alignerPath,"mem",specialsParam,"-t",str(core),args.genome, "{}"]),cmdGenerateBed,cmdGenerateBedgraph
    elif args.aligner == 'minimap2':
        if args.alignerSpecialsParam == None:
            specialsParam = "-x map-ont -a --secondary=no --MD -L"
        else:
            specialsParam = args.alignerSpecialsParam
        if args.alignerPath == None:
            alignerPath = args.aligner
        else:
            alignerPath = args.alignerPath
        return " ".join([alignerPath,args.genome,"{}",specialsParam,"-t",str(core)]),cmdGenerateBed,cmdGenerateBedgraph

def align_insert(fileSeq,outFile):
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(" ".join([baseAlignement.format(fileSeq),">",outFile]))
    tools.append_File_text(outFile+'.log', [stderr])
    if exitCode != 0:
        raise Exception("Alignement error{}: {}".format(str(exitCode),stderr))

def extract_stats_align_sam(file):
    header = True
    whereAlign = {}
    try:
        with open(file) as alignSam:
            while header:
                line = next(alignSam)
                if line[0] != '@':
                    header = False
                    cutedLine = line.split("\t")
                    whereAlign[cutedLine[2]] = 1
            while True:
                cutedLine = next(alignSam).split("\t")
                if cutedLine[2] in whereAlign:
                    whereAlign[cutedLine[2]] += 1
                else:
                    whereAlign[cutedLine[2]] = 1
    except:
        return whereAlign
    
def align_insert_barcode(outFolder,barcodeName):
    global statsAlignement
    global statsAlignementLock
    try:
        fileAlign = os.path.join(outFolder,barcodeName+".sam")
        with alignment:
            align_insert(os.path.join(outFolder,barcodeName+'_cut.fastq'),fileAlign)
            stdout, stderr, exitCode = tools.launch_cmdExt_Shell(cmdGenerateBed.replace("fileAlign",fileAlign).replace("barcodeName",barcodeName).replace("outFolder",outFolder))
            tools.test_exit_code_popen(stderr,exitCode,"align_insert_barcode during convertion sam to bed",exitScript=False)
            try:
                varian_detector(barcodeName)
            except Exception as e:
                tools.print_screen_error("{} give {}".format(barcodeName,e))
        try:
            stdout, stderr, exitCode = tools.launch_cmdExt_Shell(cmdGenerateBedgraph.replace("barcodeName",barcodeName).replace("outFolder",outFolder))
        except:
            tools.test_exit_code_popen(stderr,exitCode,"align_insert_barcode during convertion bed to bedgraph",exitScript=False)
        statsAlign = extract_stats_align_sam(fileAlign)
        with statsAlignementLock:
            if barcodeName in statsAlignement:
                for cible,count in statsAlign.items():
                    if cible in statsAlignement[barcodeName]:
                        statsAlignement[barcodeName][cible] += count
                    else:
                        statsAlignement[barcodeName][cible] = count
            else:
                statsAlignement[barcodeName] = statsAlign
    except:
        tools.print_screen_error("Error during the alignment of {}".format(barcodeName))
        return False
    
def align_insert_barcode_light(outFolder,barcodeName):
    global statsAlignement
    global statsAlignementLock
    try:
        fileAlign = os.path.join(outFolder,barcodeName+".sam")
        with alignment:
            align_insert(os.path.join(outFolder,barcodeName+'_cut.fastq'),fileAlign)
        statsAlign = extract_stats_align_sam(fileAlign)
        with statsAlignementLock:
            if barcodeName in statsAlignement:
                for cible,count in statsAlign.items():
                    if cible in statsAlignement[barcodeName]:
                        statsAlignement[barcodeName][cible] += count
                    else:
                        statsAlignement[barcodeName][cible] = count
            else:
                statsAlignement[barcodeName] = statsAlign
    except:
        return False

def align_inserts_for_visualise_in_realtime(folderRun,barcodeName):
    global statsVisualise
    global statsVisualiseLock
    global baseAlignementGraphic
    try:
        fileAlign = os.path.join(outVisualisationGraphicTmpFile,barcodeName+".sam")
        with alignment:
            stdout,stderr,exitCode = tools.launch_cmdExt_Shell(" ".join([baseAlignementGraphic.format(os.path.join(folderRun,barcodeName+'_cut.fastq')),">",fileAlign]))
        tools.append_File_text(fileAlign+'.log', [stderr])
        if exitCode != 0:
            raise Exception("Alignement error{}: {}".format(str(exitCode),stderr))
        statsAlign = extract_stats_align_sam(fileAlign)
        with statsVisualiseLock:
            if barcodeName in statsVisualise:
                for cible,count in statsAlign.items():
                    if cible in statsVisualise[barcodeName]:
                        statsVisualise[barcodeName][cible] += count
                    else:
                        statsVisualise[barcodeName][cible] = count
            else:
                statsVisualise[barcodeName] = statsAlign
    except:
        pass
    
def generate_bedgraph(outFolder,barcodeName):
    try:
        fileAlign = os.path.join(outFolder,barcodeName+".sam")
        with alignment:
            stdout, stderr, exitCode = tools.launch_cmdExt_Shell(cmdGenerateBed.replace("fileAlign",fileAlign).replace("barcodeName",barcodeName).replace("outFolder",outFolder))
            tools.test_exit_code_popen(stderr,exitCode,"generate_bedgraph during convertion sam to bed",exitScript=False)
            try:
                varian_detector(barcodeName)
            except Exception as e:
                tools.print_screen_error("{} give {}".format(barcodeName,e))
        tools.launch_cmdExt_Shell(cmdGenerateBedgraph.replace("barcodeName",barcodeName).replace("outFolder",outFolder))
    except:
        return False

def fusion_bed(bed1,bed2,out,nameOut):
    cmd = " ".join(["cat", bed1, bed2, "|", bedops, "--everything - >", os.path.join(out,nameOut+".bed")])
    stdout, stderr, exitCode = tools.launch_cmdExt_Shell(cmd)
    tools.test_exit_code_popen(stderr,exitCode,"fusion_bed",exitScript=False)
    
def merge_sam(firstSam,otherSam,out,nameOut):
    with open(os.path.join(out,nameOut+".sam"),'w') as handleOutFile:
        with open(firstSam) as handleFirstSam:
            line = handleFirstSam.readline()
            while line != '' and '@' == line[0]:
                if 'SQ' == line[1:3]:
                    handleOutFile.write(line)
                line = handleFirstSam.readline()
            while line != '':
                handleOutFile.write(line)
                line = handleFirstSam.readline()
        for samFile in otherSam:
            with open(samFile) as handleSam:
                line = handleSam.readline()
                while line != '' and '@' == line[0]:
                    line = handleSam.readline()
                while line != '':
                    handleOutFile.write(line)
                    line = handleSam.readline()

def generate_varian_detector_base_command():
    outVariant = os.path.join(outFolder,"variant")
    tools.make_dirs(outVariant)
    outVariantVisual = os.path.join(outVariant,"bamToVisualise")
    tools.make_dirs(outVariantVisual)
    #if args.vargenome != None and args.vargenome != args.genome:
    #variantCommand1 = " ".join([baseAlignement.replace(args.genome,args.vargenome).format(os.path.join(outComplex,'barcodeName_cut.fastq')),"| {} sort --output-fmt SAM -o {} -@ {} -T {}".format(samtools,os.path.join(outComplex,"barcodeName_reads_to_ref.sam"),core,os.path.join(outComplex,"barcodeName_tempPorar"))])
    #genome = args.vargenome
    #else:
    variantCommand1 = " ".join(['cat {}'.format(os.path.join(outComplex,"barcodeName.sam")),"| {} sort --output-fmt SAM -o {} -@ {} -T {}".format(samtools,os.path.join(outComplex,"barcodeName_reads_to_ref.sam"),core,os.path.join(outComplex,"barcodeName_tempPorar"))])
    genome = args.genome
    variantCommand2 = "{} --include-unpolished --no-trimming -q -1 -t {} {} {} {} > {}".format(racon,core,os.path.join(outComplex,'barcodeName_cut.fastq'),os.path.join(outComplex,"barcodeName_reads_to_ref.sam"),genome,os.path.join(outComplex,"barcodeName_racon.fa"))
    variantCommand3 = "{} sort {} --output-fmt BAM -o {} -@ {} -T {}".format(samtools,os.path.join(outComplex,"barcodeName_reads_to_ref.sam"),os.path.join(outVariantVisual,"barcodeName_reads_to_ref.bam"),core,os.path.join(outComplex,"barcodeName_tempPorar"))
    variantCommand4 = "{} index {} -@ {}".format(samtools,os.path.join(outVariantVisual,"barcodeName_reads_to_ref.bam"),core)
    variantCommand5 = "{} {} {} -x map-ont -t {} -a --secondary=no --MD -L | {} sort --output-fmt BAM -o {} -@ {} -T {}".format(minimap2,os.path.join(outComplex,"barcodeName_racon.fa"),os.path.join(outComplex,'barcodeName_cut.fastq'),core,samtools,os.path.join(outComplex,"barcodeName_reads_to_draft.bam"),core,os.path.join(outComplex,"barcodeName_tempPorar"))
    variantCommand6 = "{} index {} -@ {}".format(samtools,os.path.join(outComplex,"barcodeName_reads_to_draft.bam"),core)
    variantCommand7 = "{} consensus {} {} --model r941_min_fast_g303 --threads {}".format(medaka,os.path.join(outComplex,"barcodeName_reads_to_draft.bam"),os.path.join(outComplex,"barcodeName_consensus_probs.hdf"),core)
    variantCommand8 = "{} stitch --threads {} {} {} {}".format(medaka,core,os.path.join(outComplex,"barcodeName_consensus_probs.hdf"),os.path.join(outComplex,"barcodeName_racon.fa"),os.path.join(outComplex,"barcodeName_consensus.fasta"))
    variantCommand9 = "{} tools consensus2vcf {} {} --out_prefix {}".format(medaka,os.path.join(outComplex,"barcodeName_consensus.fasta"),genome,os.path.join(outComplex,"barcodeName_consensus_to_ref"))
    variantCommand10 = "{} tools annotate {} {} {} {} --dpsp".format(medaka,os.path.join(outComplex,"barcodeName_consensus_to_ref.vcf"),genome,os.path.join(outVariantVisual,"barcodeName_reads_to_ref.bam"),os.path.join(outComplex,"barcodeName_reads_to_ref.vcf"))
    variantCommand11 = 'awk -v mincov="'+str(args.mincov)+'''" 'BEGIN{print "#CHROM POS REF ALT COV %REF %ALT"}{c=substr($1,1,1);if(c!="#"){split($8,a,";");for(i in a){split(a[i],b,"=");if(b[1]=="DPSP"){cov=b[2]}if(b[1]=="SR"){split(b[2],d,",");s=d[1]+d[2]+d[3]+d[4];if(s>0){pref=int((d[1]+d[2])/(s)*10000)/100;palt=int((d[3]+d[4])/(s)*10000)/100}}};if(cov>mincov){print $1" "$2" "$4" "$5" "cov" "pref" "palt}}}' '''+'{} > {}'.format(os.path.join(outComplex,"barcodeName_reads_to_ref.vcf"),os.path.join(outVariant,"barcodeName_variation.txt"))
    
    return variantCommand1,variantCommand2,variantCommand3,variantCommand4,variantCommand5,variantCommand6,variantCommand7,variantCommand8,variantCommand9,variantCommand10,variantCommand11
                    
def varian_detector(barcodeName):
    '''
    Created on 15 mars 2021
    
    @author: Stéfan
    @adaptation: Francois STUDER
    
    '''
    #minimap2 $ref $fastq -x map-ont -t 1 -a --secondary=no --MD -L | samtools sort --output-fmt SAM -o ${outdir}/reads_to_ref.sam -@ 1
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand1.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 1")
    #racon --include-unpolished --no-trimming -q -1 -t 1 $fastq ${outdir}/reads_to_ref.sam $ref > ${outdir}/racon.fa
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand2.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 2")
    #samtools sort ${outdir}/reads_to_ref.sam --output-fmt BAM -o ${outdir}/reads_to_ref.bam -@ 1
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand3.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 3")
    #samtools index ${outdir}/reads_to_ref.bam -@ 1
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand4.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 4")
    #minimap2 ${outdir}/racon.fa $fastq -x map-ont -t 1 -a --secondary=no --MD -L | samtools sort --output-fmt BAM -o ${outdir}/reads_to_draft.bam -@ 1
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand5.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 5")
    #samtools index ${outdir}/reads_to_draft.bam -@ 1
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand6.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 6")
    #medaka consensus ${outdir}/reads_to_draft.bam ${outdir}/consensus_probs.hdf --model r941_min_fast_g303
    if tools.file_exists("{}/{}_consensus_probs.hdf".format(outComplex,barcodeName)):
        tools.file_remove(outComplex, barcodeName+"_consensus_probs.hdf")
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand7.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 7")
    #medaka stitch ${outdir}/consensus_probs.hdf ${outdir}/racon.fa ${outdir}/consensus.fasta
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand8.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 8")
    #medaka tools consensus2vcf ${outdir}/consensus.fasta $ref --out_prefix ${outdir}/consensus_to_ref
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand9.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 9")
    #medaka tools annotate ${outdir}/consensus_to_ref.vcf $ref ${outdir}/reads_to_ref.bam ${outdir}/reads_to_ref.vcf --dpsp
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand10.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 10",exitScript=False)
    #awk -v mincov="$mincov" 'BEGIN{print "#CHROM POS REF ALT COV %REF %ALT"}{c=substr($1,1,1);if(c!="#"){split($8,a,";");for(i in a){split(a[i],b,"=");if(b[1]=="DPSP"){cov=b[2]}if(b[1]=="SR"){split(b[2],d,",");s=d[1]+d[2]+d[3]+d[4];if(s>0){pref=int((d[1]+d[2])/(s)*10000)/100;palt=int((d[3]+d[4])/(s)*10000)/100}}};if(cov>mincov){print $1" "$2" "$4" "$5" "cov" "pref" "palt}}}' ${outdir}/reads_to_ref.vcf > ${outdir}/variation.txt
    stdout,stderr,exitCode = tools.launch_cmdExt_Shell(variantCommand11.replace("barcodeName",barcodeName))
    tools.test_exit_code_popen(stderr,exitCode,"vairan_detector step 11")
    
def trac_sequences(fastqToAnnalyse,sequencesFiles):
    idsequences = set()
    resultats = {'withoutGuib' : []}
    with open(sequencesFiles) as sequencesFilesHandle:
        for line in sequencesFilesHandle:
            idsequences.add(line.replace("\n",""))
    listJob1 = []
    for fastq in fastqToAnnalyse:
        for record in SeqIO.parse(fastq, "fastq"):
            if record.id in idsequences:
                job = calculPool.apply_async(cut_guibson, (record,querrysGuibson,lengthMiniGuibson))
                listJob1.append(job)
                
    listSecondJob = {}
    for job in listJob1:
        result = job.get()
        if len(result[1]) == 0:
            resultats['withoutGuib'].append(result[0][0].format('fastq'))
        else:
            i = result[0][0].id.split("_")[0]
            listSecondJob[i] = []
            for record in result[0]:
                job = calculPool.apply_async(test_sequence, (record,barcodes5,barcodes3,miniLength,maxLength,distanceMini,distanceMax))
                listSecondJob[i].append(job)

    
    for i,listjob in listSecondJob.items():
        resultsBarcode = {}
        results = {'multi' : {5 : [], 3 : []}, 'noGoodLength' : {-1 : [], 1 : []}, 'noBarcode' : {5 : [], 3 : []}, 'Alone 3\'' : {-1 : [], 0 : [], 1 : []}, 'barcodes':resultsBarcode}
        resultats[i] = results
        for nameBarcode5,querry in barcodes5.items():
            for nameBarcode3,querry in barcodes3.items():
                resultsBarcode[nameBarcode5+"_"+nameBarcode3] = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
        for job in listjob:
            result = job.get()
            if result[4] == 0:
                if result[2] == 0:
                    results['noBarcode'][5].append(result[0].format('fastq'))
                elif result[2] == 1:
                    if result[3] == 0:
                        results['noBarcode'][3].append(result[0].format('fastq'))
                    elif result[3] == 1:
                        resultsBarcode[result[1]][result[5]].append(result[0].format('fastq'))
                    else:
                        results['multi'][3].append(result[0].format('fastq'))
                else:
                    results['multi'][5].append(result[0].format('fastq'))
            else:
                results['noGoodLength'][result[4]].append(result[0].format('fastq'))
    
    if len(resultats['withoutGuib']) > 0:
        with open(os.path.join(outFolder,'withoutGuib'+'.txt'), "w") as file:
            for sequence in resultats['withoutGuib']:
                file.write(sequence+"\n")
    del resultats['withoutGuib']
    for seq,datas in resultats.items():
        file = open(os.path.join(outFolder,seq+'.txt'), "w")
        for key,data in datas.items():
            if key != 'barcodes':
                for subKey,subData in data.items():
                    if len(subData) > 0:
                        file.write(str(key)+"_"+str(subKey)+":\n")
                        for sequence in subData:
                            file.write(sequence)
            else:
                for barcode,groupData in data.items():
                    dataInside = False
                    for subKey,subData in groupData.items():
                        if len(subData) > 0:
                            dataInside = True
                    if dataInside:
                        file.write('#'*40+'\n'+barcode+":\n")
                        for subKey,subData in groupData.items():
                            if len(subData) > 0:
                                file.write(str(subKey)+":\n")
                                for sequence in subData:
                                    file.write(sequence)
        file.close()          

def main(fastqToAnnalyse):
    global extractedBarcode
    tools.make_dirs(outComplex)
    ParrallelFileAnnalyse = tools.init_list_thread(parrallelFile)
    for fastq in fastqToAnnalyse:
        tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=annalyse_file, args=(fastq,outComplex)))
    tools.wait_end_thread_list(ParrallelFileAnnalyse)
    
    for barcode in extractedBarcode:
        tools.auto_add_thread(ParrallelFileAnnalyse,1,Thread(target=align_insert_barcode, args=(outComplex,barcode)))

    write_finals_stats(ParrallelFileAnnalyse)
    
def get_file_from_SSH(SSHparam,setFileAlreadyDL,outFastQ,errorSSH):
    fastqToAnnalyse = []
    try:
        sshConnexion = tools.connect_SSH(SSHparam)
        sftp = sshConnexion.open_sftp()
        listFile = set(sftp.listdir(args.folder)).difference(setFileAlreadyDL)
        for fileName in listFile:
            fileToDL = os.path.join(args.folder,fileName)
            if ((time() - sftp.stat(fileToDL).st_mtime)/60) > 10:
                saveFolder = os.path.join(outFastQ,fileName)
                sftp.get(fileToDL,saveFolder)
                fastqToAnnalyse.append(saveFolder)
                setFileAlreadyDL.add(fileName)
        sftp.close()
        sshConnexion.close()
    except Exception as e:
        errorSSH += 1
        tools.print_screen_error(str(e))
    else:
        errorSSH = 0
    return fastqToAnnalyse,errorSSH

def get_file_from_system(setFileAlreadyAnnalysed,error):
    fastqToAnnalyse = []
    return fastqToAnnalyse,error

def realtime_annalyse():
    global extractedBarcode
    global begin
    
    global listStatsFiles
    global timePoints
    global chromName
    global inTimeStatsAlignement
    global colorAlignmentBC
    global colorAlignmentBCNotUsable
    global otherGenomeToVisualise
    listStatsFiles = []
    timePoints = set()
    chromName = set()
    inTimeStatsAlignement = {}
    colorAlignmentBC = {}
    colorAlignmentBCNotUsable = set()
    
    if args.visugenome != None and args.visugenome != args.genome:
        global outVisualisationGraphic
        global outVisualisationGraphicTmpFile
        global baseAlignementGraphic
        global statsVisualiseLock
        global statsVisualise
        global inTimeStatsAlignementVisualise
        global chromNameVisualise
        otherGenomeToVisualise = True
        outVisualisationGraphic = os.path.join(outFolder,'inTimeGraphic')
        outVisualisationGraphicTmpFile = os.path.join(args.tmp,name,'inTimeGraphic')
        baseAlignementGraphic = baseAlignement.replace(args.genome, args.visugenome)
        tools.make_dirs(outVisualisationGraphic)
        tools.make_dirs(outVisualisationGraphicTmpFile)
        statsVisualiseLock = RLock()
        statsVisualise = {}
        inTimeStatsAlignementVisualise = {}
        chromNameVisualise = set()
    else:
        otherGenomeToVisualise = False
    
    outTemporar = os.path.join(args.tmp,name,'temps')
    outCurentStats = os.path.join(outFolder,'inTimeStats')
    tools.make_dirs(outTemporar)
    tools.make_dirs(outCurentStats)
    ParrallelFileAnnalyse = tools.init_list_thread(parrallelFile)
    newFiles = True
    setFileAlreadyAnnalysed = set()
    runNumber = 1
    lastVariantDetection = 1
    alreadyAnnalyse = False
    timeSinceLastGet = timeBetweenCheck
    allExtractedBarcode = set()
    printStats = tools.initiate()
    mergeFiles = tools.initiate()
    #extractVariant = tools.initiate()
    printStats.start()
    mergeFiles.start()
    #extractVariant.start()
    if args.ext != None:
        SSHparam = tools.config_loader(args.ext,"SSH")
        outFastQ = os.path.join(outTemporar,'fastq')
        tools.make_dirs(outFastQ)
    error = 0
    while newFiles:
        if timeSinceLastGet < timeBetweenCheck:
            sleep((timeBetweenCheck-timeSinceLastGet)*60)
        if args.ext == None:
            fastqToAnnalyse,error = get_file_from_system(setFileAlreadyAnnalysed,error)
        else:
            fastqToAnnalyse,error = get_file_from_SSH(SSHparam,setFileAlreadyAnnalysed,outFastQ,error)
        lastGet = time()
        if len(fastqToAnnalyse) > 0:
            extractedBarcode = set()
            folderRun = os.path.join(outTemporar,str(runNumber))
            tools.make_dirs(folderRun)
            for fastq in fastqToAnnalyse:
                tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=annalyse_file, args=(fastq,folderRun)))
            tools.wait_end_thread_list(ParrallelFileAnnalyse)
    
            for barcode in extractedBarcode:
                tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=align_insert_barcode_light, args=(folderRun,barcode)))
                if otherGenomeToVisualise:
                    tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=align_inserts_for_visualise_in_realtime, args=(folderRun,barcode)))
            allExtractedBarcode.update(extractedBarcode)
            timePoint = ((time()-begin)/60/60)
            tools.wait_end_thread_list(ParrallelFileAnnalyse)
            mergeFiles.join()
            mergeFiles = Thread(target=merge_files, args=(folderRun,outComplex,extractedBarcode))
            mergeFiles.start()
            printStats.join()
            timePoints.add(timePoint)
            printStats = Thread(target=print_curent_stats, args=(outCurentStats,timePoint))
            printStats.start()
            lastProcessFile = time()
            alreadyAnnalyse = True
            if (runNumber - lastVariantDetection) > 2:
                mergeFiles.join()
                for barcode in allExtractedBarcode:
                    try:
                        varian_detector(barcode)
                    except Exception as e:
                        pass
                        #tools.print_screen_error("{} give {}".format(barcode,e))
                lastVariantDetection = runNumber
            runNumber += 1
        elif alreadyAnnalyse:
            if (time()-lastProcessFile)/60 > args.timeEndCheck:
                newFiles = False
            elif error > 20:
                newFiles = False
                
        timeSinceLastGet = (time()-lastGet)/60
    
    mergeFiles.join()
    if args.dev:
        for barcode in allExtractedBarcode:
            tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=generate_bedgraph, args=(outComplex,barcode)))
        tools.move_dir(outTemporar, os.path.join(devFolder,'temps'))
    else:
        for barcode in allExtractedBarcode:
            try:
                varian_detector(barcode)
            except Exception as e:
                pass
    extractedBarcode = allExtractedBarcode
    finalStats = Thread(target=write_finals_stats, args=(ParrallelFileAnnalyse,))
    finalStats.start()

    printStats.join()
    if (not args.dev):
        for statFile in listStatsFiles:
            tools.file_remove(outCurentStats, statFile)
    
    timePoints = list(timePoints)
    timePoints.sort()
    matrixAlign = Thread(target=save_matrix_time_point_align_stats, args=(outCurentStats,timePoints,inTimeStatsAlignement,chromName))
    matrixAlign.start()
    
    if otherGenomeToVisualise:
        save_matrix_time_point_align_stats(outVisualisationGraphic,timePoints,inTimeStatsAlignementVisualise,chromNameVisualise)
        finalStats.join()
        with open(os.path.join(outFolder,'summary.tsv'), "a") as litleStatsFile:
            litleStatsFile.write('#'*30)
            litleStatsFile.write('\nAlign statistics in visualise\n')
            write_finals_stats_align_stats(litleStatsFile,statsVisualise)
    else:
        finalStats.join()
    matrixAlign.join()
    
def save_matrix_time_point_align_stats(outMatrix,timePoints,inTimeStatsAlignementToSave,chromName):
    for chrom in chromName:
        with open(os.path.join(outMatrix,chrom+".tsv"), 'w') as chromFile:
            for timePoint in timePoints:
                chromFile.write("\t"+str(round(timePoint, 2)))
            chromFile.write('\n')
            for barcode,dataBarcode in inTimeStatsAlignementToSave.items():
                chromFile.write(barcode)
                for timePoint in timePoints:
                    chromFile.write("\t")
                    if timePoint in dataBarcode and chrom in dataBarcode[timePoint]:
                        chromFile.write(str(dataBarcode[timePoint][chrom]))
                    else:
                        chromFile.write('0')
                chromFile.write('\n')
    
def print_curent_stats(outFolder,timePoint):
    global stats
    global statsShort
    global distanceMini
    global distanceMax
    global miniLength
    global maxLength
    global errorBC
    global errorGuib
    global statsBarcodesGood
    global statsAlignement
    global otherGenomeToVisualise
    if otherGenomeToVisualise:
        global statsVisualise
    
    with writeStats:
        statsLocal = statsShort.copy()
        statsBarcodesGoodLocal = statsBarcodesGood.copy()
        statsAlignementLocal = statsAlignement.copy()
        if otherGenomeToVisualise:
            statsVisualiseLocal = statsVisualise.copy()
    global listStatsFiles
    listStatsFiles.append('summary_{}.tsv'.format(datetime.now().strftime("%Y-%m-%d_%H:%M:%S")))
    if (not args.dev):
        if len(listStatsFiles) > 3:
            tools.file_remove(outFolder, listStatsFiles.pop(0))
                
    currentStatsFile = open(os.path.join(outFolder,listStatsFiles[-1]), "w")
    TotalSeq = statsLocal['Total Seq']
    currentStatsFile.write('\tbarcode error(s) '+str(errorBC)+', guibson error(s) '+str(errorGuib)+', insert minimum length '+str(distanceMini)+', insert maximum length '+str(distanceMax)+', sequence length minimum '+str(miniLength)+', sequence max length '+str(maxLength)+'\n')
    currentStatsFile.write('Seq before cut\t'+str(statsLocal['Seq before cut'])+'\n')
    currentStatsFile.write('Amplification Factor\t'+str(TotalSeq/statsLocal['Seq before cut'])+'\n')
    currentStatsFile.write('Total Seq'+'\t'+str(TotalSeq)+'\t'+str(TotalSeq/TotalSeq*100)+'\n')
    del statsLocal['Seq before cut']
    del statsLocal['Total Seq']
    currentStatsFile.write('#'*30+'\n')
    for key,value in statsLocal.items():
        if value != 0:
            currentStatsFile.write(key+'\t'+str(value)+'\t'+str(value/TotalSeq*100)+'\n')
    currentStatsFile.write('#'*30+'\n')
    for key,value in statsBarcodesGoodLocal.items():
        if value != 0:
            currentStatsFile.write(key+'\t'+str(value)+'\t'+str(value/TotalSeq*100)+'\n')
            
    currentStatsFile.write('#'*30)
    currentStatsFile.write('\nAlign statistics\n')
    
    global chromName
    global inTimeStatsAlignement
    print_curent_stats_align_stats(currentStatsFile,timePoint,statsAlignementLocal,inTimeStatsAlignement,chromName)
    
    if otherGenomeToVisualise:
        printPlot = Thread(target=print_plot_align_stats, args=(outFolder,inTimeStatsAlignement,chromName))
        printPlot.start()
        global outVisualisationGraphic
        global outVisualisationGraphicTmpFile
        global inTimeStatsAlignementVisualise
        global chromNameVisualise
        currentStatsFile.write('#'*30)
        currentStatsFile.write('\nAlign statistics in visualise\n')
        print_curent_stats_align_stats(currentStatsFile,timePoint,statsVisualiseLocal,inTimeStatsAlignementVisualise,chromNameVisualise)
        printPlot.join()
        print_plot_align_stats(outVisualisationGraphic,inTimeStatsAlignementVisualise,chromNameVisualise)
    else:
        print_plot_align_stats(outFolder,inTimeStatsAlignement,chromName)
    
    currentStatsFile.close()
    
def print_curent_stats_align_stats(currentStatsFile,timePoint,statsAlignementLocal,inTimeStatsAlignementToUpdate,chromName):
    listBarcodeAlign = list(statsAlignementLocal.keys())
    listBarcodeAlign.sort()
    listSample = set()
    for barcode,statsAlign in statsAlignementLocal.items():
        for sample,count in statsAlign.items():
            listSample.add(sample)
            chromName.add(sample)
    listSample = list(listSample)
    listSample.sort()
    currentStatsFile.write("Sample Name")
    for sample in listSample:
        currentStatsFile.write("\t"+sample+"\tPercent"+sample)
    currentStatsFile.write("\n")
    for barcode in listBarcodeAlign:
        total = 0
        statsAlign = statsAlignementLocal[barcode]
        for sample,count in statsAlign.items():
            total += count
        currentStatsFile.write(barcode)
        if barcode in inTimeStatsAlignementToUpdate:
            dataBarcode = inTimeStatsAlignementToUpdate[barcode]
        else:
            dataBarcode = {}
            inTimeStatsAlignementToUpdate[barcode] = dataBarcode
        dataBarcode[timePoint] = {}
        for sample in listSample:
            if sample in statsAlign:
                currentStatsFile.write("\t"+str(statsAlign[sample])+"\t"+str(statsAlign[sample]/total*100))
                dataBarcode[timePoint][sample] = statsAlign[sample]
            else:
                currentStatsFile.write("\t0\t0")
        currentStatsFile.write("\n")
    
def print_plot_align_stats(outFolder,inTimeStatsAlignementUse,chromName):
    try:
        import matplotlib.pyplot as plt
        from random import randint
        import matplotlib.patheffects as mpe
        global timePoints
        global barcodes3
        global colorAlignmentBC
        global colorAlignmentBCNotUsable
        listTimePoints = list(timePoints)
        listTimePoints.sort()
        pe1 = [mpe.Stroke(linewidth=2, foreground='black'),
                mpe.Stroke(foreground='white',alpha=1),
                mpe.Normal()]
        for chrom in chromName:
            dictDataBC = {}
            for barcode,dataBarcode in inTimeStatsAlignementUse.items():
                if barcode not in colorAlignmentBC:
                    colorForBC = '#'+str(hex(randint(1118481,16777183)))[2:]
                    while colorForBC[:-1] in colorAlignmentBCNotUsable:
                        colorForBC = '#'+str(hex(randint(1118481,16777183)))[2:]
                    colorAlignmentBC[barcode] = colorForBC[:-1]+'1'
                    colorAlignmentBCNotUsable.add(colorForBC[:-1])
                listPointBarcode = []
                barcodeOnePoint = False
                for timePoint in listTimePoints:
                    if timePoint in dataBarcode and chrom in dataBarcode[timePoint]:
                        barcodeOnePoint = True
                        listPointBarcode.append(dataBarcode[timePoint][chrom])
                    else:
                        listPointBarcode.append(0)
                if barcodeOnePoint:
                    for BC in barcodes3:
                        if BC in barcode:
                            if BC in dictDataBC:
                                dictDataBC[BC].append([barcode,listPointBarcode])
                            else:
                                dictDataBC[BC] = [[barcode,listPointBarcode]]
            
            for BC,dataBC in dictDataBC.items():
                fig, ax = plt.subplots()
                for dataBCSub in dataBC:
                    ax.plot(listTimePoints,dataBCSub[1], label=dataBCSub[0], color=colorAlignmentBC[dataBCSub[0]], lw=1.5, path_effects=pe1, alpha=0.6)
                ax.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, -0.05), borderaxespad=2.5)
                ax.set_yscale('log')
                ax.set_ylabel('Number detected sequence')
                ax.set_xlabel('time (h)')
                fig.savefig(os.path.join(outFolder,BC+'_'+chrom+'.png'),dpi=400,bbox_inches='tight')
            plt.close('all')
    except Exception as e:
        tools.print_screen_error("print_plot: {}".format(e))
    
def merge_files(folderRun,outComplex,extractedBarcode):
    extractedBarcode = extractedBarcode.copy()
    for barcode in extractedBarcode:
        outSam = os.path.join(outComplex,barcode+'.sam')
        if os.path.exists(outSam):
            with open(outSam,'a') as handleOutFile:
                with open(os.path.join(folderRun,barcode+'.sam')) as handleSam:
                    line = handleSam.readline()
                    while line != '' and '@' == line[0]:
                        line = handleSam.readline()
                    while line != '':
                        handleOutFile.write(line)
                        line = handleSam.readline()
        else:
            tools.file_copy(folderRun,barcode+'.sam',outComplex)
        tools.launch_cmdExt_Shell('cat {} >> {}'.format(os.path.join(folderRun,barcode+'_cut.fastq'),os.path.join(outComplex,barcode+'_cut.fastq')))
    tools.launch_cmdExt_Shell('cat {} >> {}'.format(os.path.join(folderRun,'statOriginalSeq.tsv'),os.path.join(outComplex,'statOriginalSeq.tsv')))
    if args.dev:
        tools.launch_cmdExt_Shell('cat {} >> {}'.format(os.path.join(folderRun,'statOriginalSeqPosGuib.tsv'),os.path.join(outComplex,'statOriginalSeqPosGuib.tsv')))
    
def remove_outliers(data, thresholdStd = 3):
    """
    This method returns all values which are farther away
    than thresholdStd standard deviationa
    """
    import numpy as np
    noOutliers=[]
    mean = np.mean(data)
    std =np.std(data)
    if std == 0:
        return data
    for y in data:
        z_score= (y - mean)/std 
        if np.abs(z_score) <= thresholdStd:
            noOutliers.append(y)
    return noOutliers

def regenerate_realtime_annalyse():
    global extractedBarcode
    global begin
    
    global listStatsFiles
    global timePoints
    global chromName
    global inTimeStatsAlignement
    global colorAlignmentBC
    global colorAlignmentBCNotUsable
    global otherGenomeToVisualise
    listStatsFiles = []
    timePoints = set()
    chromName = set()
    inTimeStatsAlignement = {}
    colorAlignmentBC = {}
    colorAlignmentBCNotUsable = set()
    
    if args.visugenome != None and args.visugenome != args.genome:
        global outVisualisationGraphic
        global outVisualisationGraphicTmpFile
        global baseAlignementGraphic
        global statsVisualiseLock
        global statsVisualise
        global inTimeStatsAlignementVisualise
        global chromNameVisualise
        otherGenomeToVisualise = True
        outVisualisationGraphic = os.path.join(outFolder,'inTimeGraphic')
        outVisualisationGraphicTmpFile = os.path.join(args.tmp,name,'inTimeGraphic')
        baseAlignementGraphic = baseAlignement.replace(args.genome, args.visugenome)
        tools.make_dirs(outVisualisationGraphic)
        tools.make_dirs(outVisualisationGraphicTmpFile)
        statsVisualiseLock = RLock()
        statsVisualise = {}
        inTimeStatsAlignementVisualise = {}
        chromNameVisualise = set()
    else:
        otherGenomeToVisualise = False
    
    outTemporar = os.path.join(args.tmp,name,'temps')
    outCurentStats = os.path.join(outFolder,'inTimeStats')
    tools.make_dirs(outTemporar)
    tools.make_dirs(outCurentStats)
    ParrallelFileAnnalyse = tools.init_list_thread(parrallelFile)
    runNumber = 1
    lastVariantDetection = 1
    allExtractedBarcode = set()
    printStats = tools.initiate()
    mergeFiles = tools.initiate()
    #extractVariant = tools.initiate()
    printStats.start()
    mergeFiles.start()
    #extractVariant.start()
    fastqToAnnalyse = {}
    if args.ext != None:
        SSHparam = tools.config_loader(args.ext,"SSH")
        outFastQ = os.path.join(outTemporar,'fastq')
        tools.make_dirs(outFastQ)
        try:
            sshConnexion = tools.connect_SSH(SSHparam)
            sftp = sshConnexion.open_sftp()
            listFile = set(sftp.listdir(args.folder))
            for fileName in listFile:
                fileToDL = os.path.join(args.folder,fileName)
                saveFolder = os.path.join(outFastQ,fileName)
                if sftp.stat(fileToDL).st_mtime in fastqToAnnalyse:
                    tools.print_screen_error("File {} have the same date as {} why ?".format(fileToDL,fastqToAnnalyse[sftp.stat(fileToDL).st_mtime]))
                else:
                    fastqToAnnalyse[sftp.stat(fileToDL).st_mtime] = saveFolder
                    sftp.get(fileToDL,saveFolder)
            sftp.close()
            sshConnexion.close()
        except Exception as e:
            tools.print_screen_error(str(e))

    else:
        pass
        #fastqToAnnalyse,error = get_file_from_system(setFileAlreadyAnnalysed,error)

    if len(fastqToAnnalyse) > 0:
        listTimeForEachFile = list(fastqToAnnalyse.keys())
        listTimeForEachFile.sort()
        firstTimePoint = listTimeForEachFile[0]
        for fileTimePoint in listTimeForEachFile:
            timePoint = ((fileTimePoint-firstTimePoint)/60/60)
            extractedBarcode = set()
            folderRun = os.path.join(outTemporar,str(runNumber))
            tools.make_dirs(folderRun)
            annalyse_file(fastqToAnnalyse[fileTimePoint],folderRun)
    
            for barcode in extractedBarcode:
                tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=align_insert_barcode_light, args=(folderRun,barcode)))
                if otherGenomeToVisualise:
                    tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=align_inserts_for_visualise_in_realtime, args=(folderRun,barcode)))
            allExtractedBarcode.update(extractedBarcode)
            tools.wait_end_thread_list(ParrallelFileAnnalyse)
            mergeFiles.join()
            mergeFiles = Thread(target=merge_files, args=(folderRun,outComplex,extractedBarcode))
            mergeFiles.start()
            printStats.join()
            timePoints.add(timePoint)
            printStats = Thread(target=print_curent_stats, args=(outCurentStats,timePoint))
            printStats.start()
            if (runNumber - lastVariantDetection) > 20:
                mergeFiles.join()
                for barcode in allExtractedBarcode:
                    try:
                        varian_detector(barcode)
                    except Exception as e:
                        pass
                        #tools.print_screen_error("{} give {}".format(barcode,e))
                lastVariantDetection = runNumber
            runNumber += 1
    
        mergeFiles.join()
        if args.dev:
            for barcode in allExtractedBarcode:
                tools.auto_add_thread(ParrallelFileAnnalyse,0.5,Thread(target=generate_bedgraph, args=(outComplex,barcode)))
            tools.move_dir(outTemporar, os.path.join(devFolder,'temps'))
        else:
            for barcode in allExtractedBarcode:
                try:
                    varian_detector(barcode)
                except Exception as e:
                    pass
        extractedBarcode = allExtractedBarcode
        finalStats = Thread(target=write_finals_stats, args=(ParrallelFileAnnalyse,))
        finalStats.start()
    
        printStats.join()
        if (not args.dev):
            for statFile in listStatsFiles:
                tools.file_remove(outCurentStats, statFile)
        
        timePoints = list(timePoints)
        timePoints.sort()
        matrixAlign = Thread(target=save_matrix_time_point_align_stats, args=(outCurentStats,timePoints,inTimeStatsAlignement,chromName))
        matrixAlign.start()
        
        if otherGenomeToVisualise:
            save_matrix_time_point_align_stats(outVisualisationGraphic,timePoints,inTimeStatsAlignementVisualise,chromNameVisualise)
            finalStats.join()
            with open(os.path.join(outFolder,'summary.tsv'), "a") as litleStatsFile:
                litleStatsFile.write('#'*30)
                litleStatsFile.write('\nAlign statistics in visualise\n')
                write_finals_stats_align_stats(litleStatsFile,statsVisualise)
        else:
            finalStats.join()
        matrixAlign.join()
            
def write_finals_stats(alignThread):
    global stats
    global statsShort
    global distanceMini
    global distanceMax
    global miniLength
    global maxLength
    global errorBC
    global errorGuib
    global statsLength
    global statsAlignement
    
    TotalSeq = statsShort['Total Seq']
    generalStatsFile = open(os.path.join(outComplex,'extract_stats.tsv'), "w")
    litleStatsFile = open(os.path.join(outFolder,'summary.tsv'), "w")
    
    generalStatsFile.write('\tbarcode error(s) '+str(errorBC)+', guibson error(s) '+str(errorGuib)+', insert minimum length '+str(distanceMini)+', insert maximum length '+str(distanceMax)+', sequence length minimum '+str(miniLength)+', sequence max length '+str(maxLength)+'\n')
    generalStatsFile.write('Seq before cut\t'+str(statsShort['Seq before cut'])+'\n')
    generalStatsFile.write('Amplification Factor\t'+str(TotalSeq/statsShort['Seq before cut'])+'\n')
    litleStatsFile.write('\tbarcode error(s) '+str(errorBC)+', guibson error(s) '+str(errorGuib)+', insert minimum length '+str(distanceMini)+', insert maximum length '+str(distanceMax)+', sequence length minimum '+str(miniLength)+', sequence max length '+str(maxLength)+'\n')
    litleStatsFile.write('Seq before cut\t'+str(statsShort['Seq before cut'])+'\n')
    litleStatsFile.write('Amplification Factor\t'+str(TotalSeq/statsShort['Seq before cut'])+'\n')
    generalStatsFile.write('Total Seq'+'\t'+str(TotalSeq)+'\t'+str(TotalSeq/TotalSeq*100)+'\n')
    litleStatsFile.write('Total Seq'+'\t'+str(TotalSeq)+'\t'+str(TotalSeq/TotalSeq*100)+'\n')
    del statsShort['Seq before cut']
    del statsShort['Total Seq']
    litleStatsFile.write('#'*30+'\n')
    for key,value in statsShort.items():
        if value != 0:
            generalStatsFile.write(key+'\t'+str(value)+'\t'+str(value/TotalSeq*100)+'\n')
            litleStatsFile.write(key+'\t'+str(value)+'\t'+str(value/TotalSeq*100)+'\n')
    
    litleStatsFile.write('#'*30+'\n')
    for key,value in statsBarcodesGood.items():
        if value != 0:
            litleStatsFile.write(key+'\t'+str(value)+'\t'+str(value/TotalSeq*100)+'\n')
            
    generalStatsFile.write('#'*30)
    for key,value in stats.items():
        if value != 0:
            generalStatsFile.write(key+'\t'+str(value)+'\t'+str(value/TotalSeq*100)+'\n')
    
    for statsName,datas in statsLength.items():
        if len(datas) > 0:
            with open(os.path.join(outComplex,statsName+'.tsv'), "w") as statsFile:
                for key,champ in datas.items():
                    statsFile.write(str(key)+"\t"+str(champ)+"\n")
    
    if len(extractedBarcode) > 1:
        litleStatsFile.write('#'*30)
        litleStatsFile.write('\nAlign statistics\n')
        listExtractedBarcode = list(extractedBarcode)
        tools.wait_end_thread_list(alignThread)
        write_finals_stats_align_stats(litleStatsFile,statsAlignement)
        listBarcodeAlign = list(statsAlignement.keys())
        for barcode in listBarcodeAlign:
            total = 0
            statsAlign = statsAlignement[barcode]
            for sample,count in statsAlign.items():
                total += count
            for sample,count in statsAlign.items():
                generalStatsFile.write(barcode+"\t"+sample+"\t"+str(count)+"\t"+str(count/total*100)+"\n")
    generalStatsFile.close()
    litleStatsFile.close()
    
    if args.dev:
        fusion_bed(os.path.join(outComplex,listExtractedBarcode[0]+'.bed'),os.path.join(outComplex,listExtractedBarcode[1]+'.bed'),outComplex,'temp1')
        i = 2
        while i < len(listExtractedBarcode):
            fusion_bed(os.path.join(outComplex,'temp'+str(i-1)+'.bed'),os.path.join(outComplex,listExtractedBarcode[i]+'.bed'),outComplex,'temp'+str(i))
            tools.file_remove(outComplex,'temp'+str(i-1)+'.bed')
            i += 1
        tools.file_rename(outComplex, 'temp'+str(i-1)+'.bed', 'AllExtract.bed')
        tools.launch_cmdExt_Shell(" ".join(["sort -k 1,1",os.path.join(outComplex,"AllExtract.bed"), "|",bedtools,"genomecov -bga -split -i - -g",args.genome+'.chromsize','>',os.path.join(devFolder,"bedgraph","AllExtract.bedgraph")]))
        
        outFinalPlot = os.path.join(outFolder,'statisticsPlot')
        tools.make_dirs(outFinalPlot)
        try:
            import matplotlib.pyplot as plt
            try:
                fig, ax = plt.subplots()
                ax.plot(statsLongueurNBGuib[0], statsLongueurNBGuib[1], 'o')
                ax.set_xscale('log')
                ax.set_ylabel('Number guibson sequence detected')
                ax.set_xlabel('Length sequence')
                fig.savefig(os.path.join(outFinalPlot,'length sequence vs number guibson detected.png'),dpi=400,bbox_inches='tight')
            except Exception as e:
                err.write("Print plot error {}".format(e))
            for statsName,datas in statsLength.items():
                if len(datas) > 0:
                    length = []
                    frequence = []
                    occurence = []
                    for key,champ in datas.items():
                        length.append(int(key))
                        frequence.append(int(champ))
                        occurence.extend([int(key)]*int(champ))
                    try:
                        fig, ax = plt.subplots()
                        plt.bar(length, height=frequence)
                        ax.set_yscale('log')
                        ax.set_xscale('log')
                        ax.set_ylabel('Frequence')
                        ax.set_xlabel('Length sequence')
                        fig.savefig(os.path.join(outFinalPlot,'Histogram_'+statsName+'.png'),dpi=400,bbox_inches='tight')
                    except Exception as e:
                        err.write("Print plot error {}".format(e))
                    try:
                        fig, ax = plt.subplots()
                        plt.boxplot(occurence, showfliers=False)
                        ax.set_ylabel('Length sequence')
                        fig.savefig(os.path.join(outFinalPlot,'Boxplot_'+statsName+'.png'),dpi=400,bbox_inches='tight')
                    except Exception as e:
                        err.write("Print plot error {}".format(e))
                    try:
                        fig, ax = plt.subplots()
                        plt.violinplot(remove_outliers(occurence))
                        fig.savefig(os.path.join(outFinalPlot,'Violonplot_'+statsName+'.png'),dpi=400,bbox_inches='tight')
                    except Exception as e:
                        err.write("Print plot error {}".format(e))
        except Exception as e:
            err.write("Error during plot preparation {}".format(e))
        tools.move_dir(outComplex, os.path.join(devFolder,'debug'))
    tools.remove_dir(os.path.join(args.tmp,name))

def write_finals_stats_align_stats(litleStatsFile,statsAlignement):
    listBarcodeAlign = list(statsAlignement.keys())
    listBarcodeAlign.sort()
    listSample = set()
    for barcode,statsAlign in statsAlignement.items():
        for sample,count in statsAlign.items():
            listSample.add(sample)
    listSample = list(listSample)
    listSample.sort()
    litleStatsFile.write("Sample Name")
    for sample in listSample:
        litleStatsFile.write("\t"+sample+"\tPercent"+sample)
    litleStatsFile.write("\n")
    for barcode in listBarcodeAlign:
        total = 0
        statsAlign = statsAlignement[barcode]
        for sample,count in statsAlign.items():
            total += count
        litleStatsFile.write(barcode)
        for sample in listSample:
            if sample in statsAlign:
                litleStatsFile.write("\t"+str(statsAlign[sample])+"\t"+str(statsAlign[sample]/total*100))
            else:
                litleStatsFile.write("\t0\t0")
        litleStatsFile.write("\n")

if __name__ == '__main__':
    try:
        begin = time()
        parser = argparse.ArgumentParser(description='This script process fastq file to extract data about monomere squences', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-f","--file", metavar='file', type=str, default="", help="File(s) you want annalyse separate with commat")
        parser.add_argument("-n","--name", metavar='name', type=str, default=None, help="Name of the assay. By default, the software uses the date.")
        parser.add_argument("--folder", metavar='folder', type=str, default=None, help="Folder who contain fastq to annalyse")
        parser.add_argument("--tmp", metavar='tmp', type=str, default="/tmp", help="Folder who contain temporary files")
        parser.add_argument("-5","--cDNAbc", metavar='cDNAbc', type=str, required=True, help=" File containing the “cDNA” barcodes’ list. Format: seq-ID\tSequence")
        parser.add_argument("-3","--PCRbc", metavar='PCRbc', type=str, required=True, help="File containing the “PCR” barcodes’ list. Format: seq-ID\tSequence")
        parser.add_argument("-c","--core", metavar='core', type=int, default=1, help="Number of cores allocated to the software")
        parser.add_argument("-o","--out", metavar='outdir', type=str, default=".", help="Folder in which the output is saved")
        parser.add_argument("-G","--gibsonE", metavar='error', type=int, default=4, help="Guibson error accepted")
        parser.add_argument("-B","--barcodeE", metavar='error', type=int, default=4, help="Barcode error accepted")
        parser.add_argument("-i","--insertMinLength", metavar='insertMinLength', type=int, default=30, help="Good insert minimum length")
        parser.add_argument("-I","--insertMaxLength", metavar='insertMaxLength', type=int, default=10000, help="Good insert maximum length")
        parser.add_argument("-l","--seqMinLength", metavar='seqMinLength', type=int, default=70, help="Good sequences minimum length between guibson sequences")
        parser.add_argument("-L","--seqMaxLength", metavar='seqMaxLength', type=int, default=10000, help="Good sequences maximum length between guibson sequences")
        parser.add_argument("--noreverse3", dest='noreverse3', default=False, action='store_true', help="If you send directly the reverse complement in 3' barcodes, activate this option")
        parser.add_argument("-a","--aligner", metavar='aligner', type=str, default="bwa", help="Aligner name")
        parser.add_argument("--alignerPath", metavar='alignerPath', type=str, default=None, help="Aligner path. By default we use the name")
        parser.add_argument("--alignerSpecialsParam", metavar='alignerSpecialsParam', type=str, default=None, help="Specials parameters for the aligner")
        parser.add_argument("-g","--genome", metavar='genome', type=str, required=True, help="Aligner genome index path")
        parser.add_argument("--realtime", dest='realtime', default=False, action='store_true', help="Put the software in realtime annalyse")
        parser.add_argument("--traceseq", metavar='traceseq', type=str, default=None, help="Sequences you want track to check the method")
        parser.add_argument("--ext", metavar='ext', type=str, default=None, help="Path to the configuration for the ssh. Use this parameter if your files are not in the computer (e.g. when connecting to the Mk1c sequencer). This method are only in realtime and regeneraterealtime option are activated.")
        parser.add_argument("--timecheck", metavar='timecheck', type=int, default=10, help="In realtime analysis, you can setup the time(m) interval between each check in the folder")
        parser.add_argument("--timeEndCheck", metavar='timeEndCheck', type=int, default=40, help="In realtime analysis, set up the time(m) between the last data collection and the end of the computation due to the absence of new files to process")
        parser.add_argument("--dev", dest='dev', default=False, action='store_true', help="Keep all temporary datas")
        parser.add_argument("--regeneraterealtime", dest='regeneraterealtime', default=False, action='store_true', help="Mimic realtime analysis from all fastq files based on the information concerning the date and time of writing")
        parser.add_argument("--mincov", metavar='mincov', type=int, default=20, help="Minimum coverage to support a variation")
        parser.add_argument("-v","--visugenome", metavar='visugenome', type=str, default=None, help="Path to genome use for visualisation in realtime. This option are more if you want check interest region.")
        args = parser.parse_args()
        
        if args.file == "" and args.folder == None:
            raise Exception("You have to choose files to annalyse or a folder where the files are save.")

        if args.core < 2:
            core = 1
            calculPool = Pool(processes=core)
            listInstance = tools.init_list_thread(int(core))
        else:
            core = args.core-1
            calculPool = Pool(processes=core)
            listInstance = tools.init_list_thread(int(core*1.5))
        parrallelFile = int(core/5)+3
        numberRecordByAnnalyse = 200
        distanceMini = args.insertMinLength
        distanceMax = args.insertMaxLength
        miniLength = args.seqMinLength
        maxLength = args.seqMaxLength
        errorBC = args.barcodeE
        errorGuib = args.gibsonE
        timeBetweenCheck = args.timecheck
        if args.name != None:
            name = args.name
        else:
            name = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        outFolder = os.path.join(args.out,name)
        tools.make_dirs(outFolder)
        if args.traceseq == None:
            outComplex = os.path.join(args.tmp,name,'debug')
            tools.make_dirs(outComplex)
            if args.dev:
                devFolder = os.path.join(outFolder,'dev')
                tools.make_dirs(devFolder)
            
            baseAlignement,cmdGenerateBed,cmdGenerateBedgraph = generate_aligner_command()
            variantCommand1,variantCommand2,variantCommand3,variantCommand4,variantCommand5,variantCommand6,variantCommand7,variantCommand8,variantCommand9,variantCommand10,variantCommand11 = generate_varian_detector_base_command()

        querrysGuibson['+'] = []
        querrysGuibson['-'] = []
        lengthMiniGuibson = len(sequenceGuibson)+args.gibsonE
        for i in range(0,args.gibsonE+1):
            querrysGuibson['+'].append(['('+sequenceGuibson+'){e<='+str(i)+'}',[]])# (?e) est ce vraiment utile lors de ce parcours ? Notre recherche tant a remplacer.
            querrysGuibson['-'].append(['('+str(Seq(sequenceGuibson).reverse_complement())+'){e<='+str(i)+'}',[]])
        with open(args.cDNAbc, "r") as bc5File:
            for datas in bc5File:
                cutDatas = datas[:-1].split("\t")
                barcodes5[cutDatas[0]] = {'+' : '(?e)('+cutDatas[1]+'){e<='+str(args.barcodeE)+'}', '-' : '(?e)('+str(Seq(cutDatas[1]).reverse_complement())+'){e<='+str(args.barcodeE)+'}'}
        with open(args.PCRbc, "r") as bc3File:
            for datas in bc3File:
                cutDatas = datas[:-1].split("\t")
                if args.noreverse3:
                    barcodes3[cutDatas[0]] = {'+' : '(?e)('+cutDatas[1]+'){e<='+str(args.barcodeE)+'}', '-' : '(?e)('+str(Seq(cutDatas[1]).reverse_complement())+'){e<='+str(args.barcodeE)+'}'}
                else:
                    barcodes3[cutDatas[0]] = {'+' : '(?e)('+str(Seq(cutDatas[1]).reverse_complement())+'){e<='+str(args.barcodeE)+'}', '-' : '(?e)('+cutDatas[1]+'){e<='+str(args.barcodeE)+'}'}
        for nameBarcode5,querry in barcodes5.items():
            for nameBarcode3,querry in barcodes3.items():
                nameBarcode = nameBarcode5+"_"+nameBarcode3
                stats['Number good unknow '+nameBarcode] = 0
                stats['Number good '+nameBarcode] = 0
                stats['Number good too close '+nameBarcode] = 0
                stats['Number good too far '+nameBarcode] = 0
                stats['Number good bad strand '+nameBarcode] = 0
                stats['Number good bad position '+nameBarcode] = 0
                statsBarcodesGood['Number good '+nameBarcode] = 0
        fastqToAnnalyse = tools.remove_empty_string_list(tools.remove_None_in_list(args.file.split(",")))
        if args.regeneraterealtime:
            regenerate_realtime_annalyse()
        elif args.realtime:
            realtime_annalyse()
        else:
            if args.folder != None:
                for file in os.listdir(args.folder):
                    if file.endswith(".fastq"):
                        fastqToAnnalyse.append(os.path.join(args.folder,file))
            if args.traceseq == None:
                main(fastqToAnnalyse)
            else:
                trac_sequences(fastqToAnnalyse,args.traceseq)
        tools.print_screen("The software take {} to treat all data.".format(str(time()-begin)))
    except Exception as e:
        err.write(str(e)+"\n")
        exit(1)
    exit(0)
