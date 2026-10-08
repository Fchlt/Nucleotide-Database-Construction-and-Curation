##### 1- GetGenefromGenome ##### 
library(seqinr)

# Path to genome directories
PATH='ncbi_dataset/ncbi_dataset/data/'
#load the genome info
genome.info=read.delim(paste(PATH,'data_summary.tsv',sep=''))
# Get a list of genome accessions (directories to search)
genomes=list.files(PATH)
genomes=genomes[-grep('json',genomes)]
genomes=genomes[-grep('tsv',genomes)]
# Load a list of search terms
Target='amoA'
Use.Alternative.Queries=T
Queries=read.delim('amoA_search_terms.txt',header = F)
strict=F

Summary.search=NULL

for (i in 1:length(genomes)){
  
  files=list.files(paste(PATH,genomes[i],sep=''))
  
  if("cds_from_genomic.fna"%in%files){
    #find the organism name
    if('Organism.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[i],'Organism.Name']}
    if('Organism.Scientific.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[i],'Organism.Scientific.Name']}
    genome.seq=seqinr::read.fasta(paste(paste(PATH,genomes[i],sep=''),'cds_from_genomic.fna',sep='/'),whole.header = T,forceDNAtolower = F,as.string = T)
    
    #if the organism name is not found, check if it's due to issues with genome names as annotated by NCBI RefSeq and GenBank
    if(identical(Organism.Name,character(0))){
      Start=unlist(strsplit(genomes[i],split='_'))[1]
      if(Start=='GCF'){genomes[i]=gsub('GCF','GCA',genomes[i])}
      if(Start=='GCA'){genomes[i]=gsub('GCA','GCF',genomes[i])}
      if('Organism.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[i],'Organism.Name']}
      if('Organism.Scientific.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[i],'Organism.Scientific.Name']}
    }
    
    GENE=NULL
    #find the target gene
    if(strict==T){GENE=grep(paste('gene=',Target,']',sep=''),names(genome.seq))}else{GENE=grep(paste('gene=',Target,sep=''),names(genome.seq))}
    #if th egene name hasn't been found, use the alternative search queries
    if(identical(GENE,integer(0))&Use.Alternative.Queries==T){
      hit=F
      index=1
      while(hit==F&index<=nrow(Queries)){
        GENE=grep(Queries[index,1],names(genome.seq))
        if(identical(GENE,integer(0))==F){hit=T}else{index=index+1}
      }
    }
    
    #if there's a hit, print the gene in th efasta file
    if (identical(GENE,integer(0))==F){
      GENE.seq=genome.seq[GENE]
      
      New.names=NULL
      Protein.names=NULL
      for(j in 1:length(GENE)){
        New.names[j]=gsub('-','_',gsub('\\)','',gsub('\\(','',gsub(' ','_',paste(Organism.Name,gsub('\\|','',gsub('lcl','',unlist(strsplit(names(GENE.seq),split = ' '))[1])))))))
        Protein.names[j]=gsub(' ','_',(gsub('-','_',unlist(strsplit(unlist(strsplit(names(GENE.seq),split = 'protein='))[2],split='] '))[1])))
        New.names[j]=paste(genomes[i],New.names[j],Protein.names[j],sep='_')
      }
      
      write.fasta(sequences = GENE.seq, names = New.names,nbchar = max(nchar(GENE.seq)), file.out=paste(Target,'.fasta',sep=''),open="a")
      
      tmp=data.frame(
        Sequence=New.names,
        Genome=genomes[i],
        Species=Organism.Name,
        Gene.found='Yes',
        N.hits=length(GENE)
      )
      
    }else{
      tmp=data.frame(
        Sequence='NA',
        Genome=genomes[i],
        Species=Organism.Name,
        Gene.found=paste(Target,'not found in genome',sep=''),
        N.hits=0)
    }
    
  }else{
    tmp=data.frame(
      Sequence='NA',
      Genome=genomes[i],
      Species=Organism.Name,
      Gene.found='No genome found',
      N.hits=0)
  }
  
  if(i==1){Summary.search=tmp}else{Summary.search=rbind(Summary.search,tmp)}
  
  print(paste('processing',(round(i/length(genomes),5)*100),'%'))
  if(i==length(genomes)){write.csv(Summary.search,paste(Target,'Summary search.csv'),row.names = F)}
}

##### 2- ORFchecker ##### 
library(Biostrings)

Target='rpoB'
sequences=readDNAStringSet('../Construction/rpoB.fasta')
Length_threshold=1999
START.list=c('ATG','GTG','TTG','CTG','ATT','ATC')

for (i in 1:length(sequences)){
  if(length(sequences[[i]])>Length_threshold){
    
    START=paste(sequences[[i]][1:3],collapse='')
    PROT=seqinr::translate(unlist(strsplit(as.character(sequences[i]),split='')))
    RC=Biostrings::reverseComplement(sequences[i])
    START.RC=paste(RC[[1]][1:3],collapse='')
    PROT.RC=seqinr::translate(unlist(strsplit(as.character(RC),split='')))
    
    ORF=ifelse('*'%in%PROT[1:(length(PROT)-1)],F,T)
    nonAmbiguous=ifelse('X'%in%PROT,F,T)
    
    ORF.RC=ifelse('*'%in%PROT.RC[1:(length(PROT.RC)-1)],F,T)
    nonAmbiguous.RC=ifelse('X'%in%PROT.RC,F,T)
    
    if(START%in%START.list & PROT[length(PROT)]=='*' & ORF & nonAmbiguous){
      writeXStringSet(sequences[i],filepath=paste(Target,paste(database,'_full_ORF.fasta',sep=''),sep='_'), append = T)
    }else if(START.RC%in%START.list & PROT.RC[length(PROT.RC)]=='*' & ORF.RC & nonAmbiguous.RC){
      writeXStringSet(RC,filepath=paste(Target,paste(database,'_full_ORF_ReverseComplement.fasta',sep=''),sep='_'), append = T)
    }else{
      writeXStringSet(sequences[i],filepath=paste(Target,paste(database,'_partial.fasta',sep=''),sep='_'), append = T)
    }
    
  }else{writeXStringSet(sequences[i],filepath=paste(Target,paste(database,'tooShort.fasta',sep=''),sep='_'), append = T)}
  
  print(paste('processing',(round(i/length(sequences),5)*100),'%'))
}

##### 3- FindSpecies ##### 
library(seqinr)

Target='rpoB'
database='NCBI_datasets'
sequences=seqinr::read.fasta(paste(Target,paste(database,'_full_ORF.fasta',sep=''),sep='_'),whole.header=T, as.string = T)
genome.list=read.delim('genomelist.txt')

Species=NULL
for (i in 1:length(sequences)){
  if(database=='ENA'){
    
    if(unlist(strsplit(names(sequences)[i],split=' '))[2]=='Candidatus'|unlist(strsplit(names(sequences)[i],split=' '))[2]=='candidatus'){
      gen=paste(unlist(strsplit(names(sequences)[i],split=' '))[2],unlist(strsplit(names(sequences)[i],split=' '))[3])
      spe=paste(unlist(strsplit(names(sequences)[i],split=' '))[2],unlist(strsplit(names(sequences)[i],split=' '))[3],unlist(strsplit(names(sequences)[i],split=' '))[4])
    }else{
      gen=unlist(strsplit(names(sequences)[i],split=' '))[2]
      spe=paste(unlist(strsplit(names(sequences)[i],split=' '))[2],unlist(strsplit(names(sequences)[i],split=' '))[3])
    }
    
  }
  
  if(database=='NCBI_nucleotide'){
    if(unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[2]=='Candidatus'|unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[2]=='candidatus'){
      gen=paste(unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[2],unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[3])
      spe=paste(unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[2],unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[3],unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[4])
    }
    gen=unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[2]
    spe=paste(unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[2],unlist(strsplit(gsub(' UNVERIFIED:','',names(sequences)[i]),split=' '))[3])
  }
  
  if(database=='MGDB'){
    abbr=unlist(strsplit(names(sequences)[i],split=':'))[1]
    gen=unlist(strsplit(genome.list[genome.list$abbrev==abbr,'organism'],split=' '))[1]
    spe=paste(unlist(strsplit(genome.list[genome.list$abbrev==abbr,'organism'],split=' '))[1],unlist(strsplit(genome.list[genome.list$abbrev==abbr,'organism'],split=' '))[2])
    
    if(is.null(gen)){gen='Unknown'}
    if(identical(spe,character(0))){spe='Unknown'}
  }
  
  if(database=='NCBI_datasets'){
    
    UNL=unlist(strsplit(names(sequences[i]),split='_'))
    if(UNL[3]=='Candidatus'){
      gen=paste(UNL[c(3,4)],collapse=' ')
      spe=paste(UNL[c(3,4,5)],collapse=' ')
    }else if('Phytoplasma'%in%UNL[3:length(UNL)]|'phytoplasma'%in%UNL[3:length(UNL)]){
      gen='Phytoplasma'
      spe='Phytoplasma sp.'
    }else{
      gen=paste(UNL[c(3)],collapse=' ')
      spe=paste(UNL[c(3,4)],collapse=' ')
    }
  }
  
  tmp=data.frame(
    seq.ID=names(sequences)[i],
    genus=gen,
    species=spe
  )
  
  if(i==1){Species=tmp}else{Species=rbind(Species,tmp)}
  
  print(paste('progress : ',i,'/',length(sequences)))
  
}
write.csv(Species, paste(Target,paste(database,'_full_ORF_species_list.csv',sep=''),sep='_') ,row.names = F)

##### 4- GenerateTax ##### 
library(seqinr)
library(taxize)
library(openxlsx)
library(Biostrings)

Target='rpoB'
database='NCBI_datasets'
species=read.csv(paste(Target,database,'full_ORF_species_list.csv',sep='_'))
Sys.setenv(ENTREZ_KEY = "3799b1a2d846461c496f0b2475f33244bf08")

NCBI.Tax=NULL
NBC=NULL
BLCA=NULL

for(i in 1:nrow(species)){
  
  GENUS=NULL
  SPECIES=NULL
  
  GENUS=species$genus[i]
  if(identical(grep('uncultured',GENUS),integer(0))==F){GENUS=gsub('uncultured ','',GENUS)}
  SPECIES=species$species[i]
  
  SPECIES=gsub(']','',(gsub('[','',SPECIES,fixed = T)),fixed = T)
  GENUS=gsub(']','',(gsub('[','',GENUS,fixed = T)),fixed = T)
  
  if(identical(grep(' NBC ',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' NBC '))[1]}
  if(identical(grep(' NBRC ',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' NBRC '))[1]}
  if(identical(grep(' NBT',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' NBT'))[1]}
  if(identical(grep(' ANB',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' ANB'))[1]}
  if(identical(grep(' NBA',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' NBA'))[1]}
  if(identical(grep(' NH8B',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' NH8B'))[1]}
  if(identical(grep(' NHP',SPECIES),integer(0))==F){SPECIES=unlist(strsplit(SPECIES,split=' NHP'))[1]}
  
  if(i==1){NCBI.Tax.backup=NULL}
  
  if(i>1&SPECIES%in%NCBI.Tax.backup$Species){
    tmp=data.frame(
      seq.ID=species$seq.ID[i],
      Domain=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'Domain'])[1],
      Phylum=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'Phylum'])[1],
      Class=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'Class'])[1],
      Order=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'Order'])[1],
      Family=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'Family'])[1],
      Genus=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'Genus'])[1],
      Species=SPECIES,
      NCBI.ID=unique(NCBI.Tax.backup[NCBI.Tax.backup$Species==SPECIES,'NCBI.ID'])[1]
    )
  }else{
    
    search=classification(SPECIES, db = 'ncbi')
    result=as.data.frame(search[1])
    
    if(nrow(result)>1){
      if('domain'%in%result[,2]){DO=result[result[,2]=='domain',1]}else{DO='Unknown'}
      if('phylum'%in%result[,2]){PH=result[result[,2]=='phylum',1]}else{PH='Unknown'}
      if('class'%in%result[,2]){CL=result[result[,2]=='class',1]}else{CL='Unknown'}
      if('order'%in%result[,2]){OR=result[result[,2]=='order',1]}else{OR='Unknown'}
      if('family'%in%result[,2]){FA=result[result[,2]=='family',1]
      }else if('order'%in%result[,2]&'genus'%in%result[,2]){
        if(grep('genus',result[,2])==grep('order',result[,2])+2){
          FA=result[median(c(grep('genus',result[,2]),grep('order',result[,2]))),1]
        }
      }else{FA='Unknown'}
      if('genus'%in%result[,2]){GE=result[result[,2]=='genus',1]}
      if('species'%in%result[,2]){SPECIES=result[result[,2]=='species',1]}
      ID=result[nrow(result),3]
    }else{
      DO='Unknown'
      PH='Unknown'
      CL='Unknown'
      OR='Unknown'
      FA='Unknown'
      GE='Unknown'
      ID='Unknown'
    }
    tmp=data.frame(
      seq.ID=species$seq.ID[i],
      Domain=DO,
      Phylum=PH,
      Class=CL,
      Order=OR,
      Family=FA,
      Genus=GENUS,
      Species=SPECIES,
      NCBI.ID=ID)
    
    #If nothing has been found, try changing the names of the unnassigned species
    if(DO=='Unknown'&identical(grep('sp.',SPECIES,fixed = T),integer(0))==F){
      
      new.spe.name=paste('uncultured',gsub(' sp.','',SPECIES,fixed=T),sep=' ')
      
      search=classification(new.spe.name, db = 'ncbi')
      result=as.data.frame(search[1])
      
      if(nrow(result)>1){
        if('domain'%in%result[,2]){DO=result[result[,2]=='domain',1]}else{DO='Unknown'}
        if('phylum'%in%result[,2]){PH=result[result[,2]=='phylum',1]}else{PH='Unknown'}
        if('class'%in%result[,2]){CL=result[result[,2]=='class',1]}else{CL='Unknown'}
        if('order'%in%result[,2]){OR=result[result[,2]=='order',1]}else{OR='Unknown'}
        if('family'%in%result[,2]){FA=result[result[,2]=='family',1]}else{FA='Unknown'}
        if('genus'%in%result[,2]){GE=result[result[,2]=='genus',1]}else{GE='Unknown'}
        ID=result[nrow(result),3]
      }else{
        new.spe.name=paste('unclassified',gsub(' sp.','',SPECIES,fixed=T),sep=' ')
        search=classification(new.spe.name, db = 'ncbi')
        result=as.data.frame(search[1])
        if(nrow(result)>1){
          if('domain'%in%result[,2]){DO=result[result[,2]=='domain',1]}else{DO='Unknown'}
          if('phylum'%in%result[,2]){PH=result[result[,2]=='phylum',1]}else{PH='Unknown'}
          if('class'%in%result[,2]){CL=result[result[,2]=='class',1]}else{CL='Unknown'}
          if('order'%in%result[,2]){OR=result[result[,2]=='order',1]}else{OR='Unknown'}
          if('family'%in%result[,2]){FA=result[result[,2]=='family',1]}else{FA='Unknown'}
          if('genus'%in%result[,2]){GE=result[result[,2]=='genus',1]}else{GE='Unknown'}
          ID=result[nrow(result),3]
        }else{
          DO='Unknown'
          PH='Unknown'
          CL='Unknown'
          OR='Unknown'
          FA='Unknown'
          GE='Unknown'
          ID='Unknown'
        }
      }
      tmp=data.frame(
        seq.ID=species$seq.ID[i],
        Domain=DO,
        Phylum=PH,
        Class=CL,
        Order=OR,
        Family=FA,
        Genus=GENUS,
        Species=SPECIES,
        NCBI.ID=ID)
    }
    
    #Finally, if nothing has worked, search for genus name
    if(DO=='Unknown'){
      search=classification(GENUS, db = 'ncbi')
      result=as.data.frame(search[1])
      if(nrow(result)>1){
        if('domain'%in%result[,2]){DO=result[result[,2]=='domain',1]}else{DO='Unknown'}
        if('phylum'%in%result[,2]){PH=result[result[,2]=='phylum',1]}else{PH='Unknown'}
        if('class'%in%result[,2]){CL=result[result[,2]=='class',1]}else{CL='Unknown'}
        if('order'%in%result[,2]){OR=result[result[,2]=='order',1]}else{OR='Unknown'}
        if('family'%in%result[,2]){FA=result[result[,2]=='family',1]
        }else if('order'%in%result[,2]&'genus'%in%result[,2]){
          if(grep('genus',result[,2])==grep('order',result[,2])+2){
            FA=result[median(c(grep('genus',result[,2]),grep('order',result[,2]))),1]
          }
        }else{FA='Unknown'}
        if('genus'%in%result[,2]){GE=result[result[,2]=='genus',1]}
        ID=result[nrow(result),3]
      }else{
        DO='Unknown'
        PH='Unknown'
        CL='Unknown'
        OR='Unknown'
        FA='Unknown'
        GE='Unknown'
        ID='Unknown'
      }
      tmp=data.frame(
        seq.ID=species$seq.ID[i],
        Domain=DO,
        Phylum=PH,
        Class=CL,
        Order=OR,
        Family=FA,
        Genus=GENUS,
        Species=SPECIES,
        NCBI.ID=ID)
    }
    
  }
  
  tmp.NBC=data.frame(
    seq.ID=species$seq.ID[i],
    Tax.NBC=paste(
      paste('k__',DO,sep=''),
      paste('p__',PH,sep=''),
      paste('c__',CL,sep=''),
      paste('o__',OR,sep=''),
      paste('f__',FA,sep=''),
      paste('g__',GENUS,sep=''),
      paste('s__',SPECIES,sep=''),sep='; '))
  
  tmp.BLCA=data.frame(
    seq.ID=species$seq.ID[i],
    Tax.BLCA=paste(
      paste('species:',DO,sep=''),
      paste('genus:',PH,sep=''),
      paste('family:',CL,sep=''),
      paste('order:',OR,sep=''),
      paste('class:',FA,sep=''),
      paste('phylum:',GENUS,sep=''),
      paste('superkingdom:',SPECIES,sep=''),sep=';')) 
  
  if(i==1){NCBI.Tax=tmp;NBC=tmp.NBC;BLCA=tmp.BLCA}else{NCBI.Tax=rbind(NCBI.Tax,tmp);NBC=rbind(NBC,tmp.NBC);BLCA=rbind(BLCA,tmp.BLCA)}
  
  NCBI.Tax.backup=NCBI.Tax
  
  print(paste(i,'/',nrow(species),sep=' '))
  if(i==nrow(species)){
    write.csv(NCBI.Tax,paste(Target,database,'full_ORF_Taxonomy.csv',sep='_'),row.names = F)
    write.table(NBC,paste(Target,database,'full_ORF_Taxonomy_NBC.txt',sep='_'),col.names = F,row.names = F,quote = F,sep='\t')
    write.table(BLCA,paste(Target,database,'full_ORF_Taxonomy_BLCA.txt',sep='_'),,col.names = F,row.names = F,quote = F,sep='\t')
  }
}



