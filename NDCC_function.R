NDCC=function(PATH,
              annotation='ncbi',
              Target,
              Length_threshold=0,
              Use.Alternative.Queries=F,
              Queries=NULL,
              strict=F,
              START.list=c('ATG','GTG','TTG','CTG','ATT','ATC')){
  
  require(seqinr)
  require(Biostrings)
  require(taxize)
  require(openxlsx)
  require(dbplyr)
  
  ############ 1- GetGenefromGenome ############ 
  database='NCBI_datasets'
  genomes=list.files(PATH)
  INFO=genomes[grep('.tsv',genomes)]
  genome.info=read.delim(paste(PATH,INFO,sep='/'))
  
  if(identical(grep('json',genomes),integer(0))==F){genomes=genomes[-grep('json',genomes)]}
  if(identical(grep('tsv',genomes),integer(0))==F){genomes=genomes[-grep('tsv',genomes)]}
  if(identical(grep('fna',genomes),integer(0))==F){genomes=genomes[-grep('fna',genomes)]}
  
  # find a list of unique accession 
  if(identical(grep('GCF_',genomes),integer(0))==F){
    GCF=genomes[grep('GCF_',genomes)]
    GCF_to_GCA=gsub('GCF_','GCA_',GCF)
    GCF_to_GCA=GCF_to_GCA[!GCF_to_GCA%in%genomes]
    if(length(GCF_to_GCA)>0){
      GCA_to_GCF=gsub('GCA_','GCF_',GCF_to_GCA)
      genomes=unique(c(genomes[-grep('GCF_',genomes)],GCA_to_GCF))
    }
  }
  
  print(paste('Extracting',Target,'genes from',length(genomes),'genomes'))
  
  # Load a list of search terms
  Summary.search=NULL
  for (g in 1:length(genomes)){
    files=list.files(paste(PATH,genomes[g],sep='/'))
    
    if(annotation=='ncbi'){FILE=ifelse("cds_from_genomic.fna"%in%files,T,F)}
    if(annotation=='prokka'){FILE=ifelse(identical(grep('.ffn',files,fixed=T),integer(0)),F,T)}
    
    if(FILE){
      #find the organism name
      if('Organism.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[g],'Organism.Name']}
      if('Organism.Scientific.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[g],'Organism.Scientific.Name']}
      if(annotation=='ncbi'){genome.seq=seqinr::read.fasta(paste(paste(PATH,genomes[g],sep=''),'cds_from_genomic.fna',sep='/'),whole.header = T,forceDNAtolower = F,as.string = T)}
      if(annotation=='prokka'){genome.seq=seqinr::read.fasta(paste(paste(PATH,genomes[g],sep=''),files[grep('ffn',files)],sep='/'),whole.header = T,forceDNAtolower = F,as.string = T)}
      
      #if the organism name is not found, check if it's due to issues with genome names as annotated by NCBI RefSeq and GenBank
      if(identical(Organism.Name,character(0))){
        Start=unlist(strsplit(genomes[g],split='_'))[1]
        if(Start=='GCF'){genomes[g]=gsub('GCF','GCA',genomes[g])}
        if(Start=='GCA'){genomes[g]=gsub('GCA','GCF',genomes[g])}
        if('Organism.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[g],'Organism.Name']}
        if('Organism.Scientific.Name'%in%colnames(genome.info)){Organism.Name=genome.info[genome.info$Assembly.Accession==genomes[g],'Organism.Scientific.Name']}
      }
      
      Organism.Name=unique(Organism.Name)
      
      GENE=NULL
      #find the target gene
      if(strict==T){GENE=grep(paste('gene=',Target,']',sep=''),names(genome.seq))}else{GENE=grep(paste('gene=',Target,sep=''),names(genome.seq))}
      #if the gene name hasn't been found, use the alternative search queries
      if(identical(GENE,integer(0))&Use.Alternative.Queries==T){
        hit=F
        index=1
        while(hit==F&index<=nrow(Queries)){
          if(strict==T){
            names(genome.seq)=paste(names(genome.seq),']',sep='')
            GENE=grep(paste(Queries[index,1],']',sep=''),names(genome.seq))
          }else{
            GENE=grep(Queries[index,1],names(genome.seq))
          }
          if(identical(GENE,integer(0))==F){hit=T}else{index=index+1}
        }
      }
      
      #if there's a hit, print the gene in th fasta file
      if (identical(GENE,integer(0))==F){
        GENE.seq=genome.seq[GENE]
        
        New.names=NULL
        Protein.names=NULL
        for(j in 1:length(GENE)){
          New.names[j]=gsub(']','',gsub('-','_',gsub('\\)','',gsub('\\(','',gsub(' ','_',paste(Organism.Name,gsub('\\|','',gsub('lcl','',unlist(strsplit(names(GENE.seq),split = ' '))[1]))))))))
          if(annotation=='ncbi'){
            Protein.names[j]=gsub(' ','_',(gsub('-','_',unlist(strsplit(unlist(strsplit(names(GENE.seq),split = 'protein='))[2],split='] '))[1])))
            New.names[j]=paste(genomes[g],New.names[j],Protein.names[j],sep='_')
          }else{
            New.names[j]=paste(genomes[g],New.names[j],sep='_')
          }
          New.names[j]=paste('seq',paste(g,j,sep='.'),New.names[j],sep='_')
        }
        
        
        seqinr::write.fasta(sequences = GENE.seq, names = New.names,nbchar = max(nchar(GENE.seq)), file.out=paste(Target,'.fasta',sep=''),open="a")
        
        tmp=data.frame(
          Sequence=New.names,
          Genome=genomes[g],
          Species=Organism.Name,
          Gene.found='Yes',
          N.hits=length(GENE)
        )
        
      }else{
        tmp=data.frame(
          Sequence='NA',
          Genome=genomes[g],
          Species=Organism.Name,
          Gene.found=paste(Target,' not found in genome',sep=''),
          N.hits=0)
      }
      
    }else{
      tmp=data.frame(
        Sequence='NA',
        Genome=genomes[g],
        Species=Organism.Name,
        Gene.found='No genome found',
        N.hits=0)
    }
    
    if(g==1){Summary.search=tmp}else{Summary.search=rbind(Summary.search,tmp)}
    
  }
  
  Summary.search.sub=dplyr::distinct(Summary.search[,c('Genome','Species','Gene.found','N.hits')], Genome, .keep_all = T)
  write.csv(Summary.search.sub,paste(Target,'Summary search.csv'),row.names = F)
  print(paste(sum(Summary.search.sub$N.hits),Target,'genes extracted from',length(genomes),'genomes'))
  if(nrow(Summary.search.sub[Summary.search.sub$N.hits==0,])>0){print(paste('No',Target,'gene found in',nrow(Summary.search.sub[Summary.search.sub$N.hits==0,]),'genomes'))}
  
  #############################################
  
  ############ 2- ORFchecker ############
  print('Extracting full-length genes')
  
  sequences=Biostrings::readDNAStringSet(paste(Target,'.fasta',sep=''))
  for (s in 1:length(sequences)){
    
    START=paste(sequences[[s]][1:3],collapse='')
    PROT=seqinr::translate(unlist(strsplit(as.character(sequences[s]),split='')))
    RC=Biostrings::reverseComplement(sequences[s])
    START.RC=paste(RC[[1]][1:3],collapse='')
    PROT.RC=seqinr::translate(unlist(strsplit(as.character(RC),split='')))
    
    ORF=ifelse('*'%in%PROT[1:(length(PROT)-1)],F,T)
    nonAmbiguous=ifelse('X'%in%PROT,F,T)
    
    ORF.RC=ifelse('*'%in%PROT.RC[1:(length(PROT.RC)-1)],F,T)
    nonAmbiguous.RC=ifelse('X'%in%PROT.RC,F,T)
    
    if(START%in%START.list & PROT[length(PROT)]=='*' & ORF & nonAmbiguous & length(sequences[[s]])>Length_threshold){
      Biostrings::writeXStringSet(sequences[s],filepath=paste(Target,paste(database,'_full_ORF.fasta',sep=''),sep='_'), append = T)
    }else if(START.RC%in%START.list & PROT.RC[length(PROT.RC)]=='*' & ORF.RC & nonAmbiguous.RC & length(sequences[[s]])>Length_threshold){
      Biostrings::writeXStringSet(RC,filepath=paste(Target,paste(database,'_full_ORF.fasta',sep=''),sep='_'), append = T)
    }else if(ORF==F){
      Biostrings::writeXStringSet(sequences[s],filepath=paste(Target,paste(database,'_NON_ORF.fasta',sep=''),sep='_'), append = T)
    }else if(ORF){
      Biostrings::writeXStringSet(sequences[s],filepath=paste(Target,paste(database,'_partial_ORF.fasta',sep=''),sep='_'), append = T)
    }else if(ORF.RC){
      Biostrings::writeXStringSet(RC,filepath=paste(Target,paste(database,'_partial_ORF.fasta',sep=''),sep='_'), append = T)
    }
  }
  #check that full-length sequences have been found
  files.2=list.files()
  if(!paste(Target,paste(database,'_full_ORF.fasta',sep=''),sep='_')%in%files.2){
    stop(paste('No full ORF',Target,'found'))
  }
  
  sequences.f=seqinr::read.fasta(paste(Target,paste(database,'_full_ORF.fasta',sep=''),sep='_'),whole.header=T, as.string = T)
  print(paste(length(sequences.f),'full ORF found out of',length(sequences),'input sequences'))
  #############################################
  
  ##### 3- FindSpecies ##### 
  Species=NULL
  for (f in 1:length(sequences.f)){
    
    UNL=unlist(strsplit(names(sequences.f)[f],split='_'))
    if(UNL[5]=='Candidatus'){
      gen=paste(UNL[c(5,6)],collapse=' ')
      spe=paste(UNL[c(5,6,7)],collapse=' ')
    }else if('Phytoplasma'%in%UNL[5:length(UNL)]|'phytoplasma'%in%UNL[5:length(UNL)]){
      gen='Phytoplasma'
      spe='Phytoplasma sp.'
    }else{
      gen=paste(UNL[c(5)],collapse=' ')
      spe=paste(UNL[c(5,6)],collapse=' ')
    }
    
    tmp=data.frame(
      seq.ID=names(sequences.f)[f],
      genus=gen,
      species=spe
    )
    
    if(f==1){Species=tmp}else{Species=rbind(Species,tmp)}
    
  }
  write.csv(Species, paste(Target,paste(database,'_full_ORF_species_list.csv',sep=''),sep='_') ,row.names = F)
  #############################################
  
  ##### 4- GenerateTax ##### 
  print('Generating reference taxonomy files')
  species=read.csv(paste(Target,paste(database,'_full_ORF_species_list.csv',sep=''),sep='_'))
  NCBI.Tax=NULL
  NBC=NULL
  BLCA=NULL
  
  for(t in 1:nrow(species)){
    
    GENUS=NULL
    SPECIES=NULL
    
    GENUS=species$genus[t]
    SPECIES=species$species[t]
    if(identical(grep('uncultured',GENUS),integer(0))==F){
      GENUS=unlist(strsplit(SPECIES,split=' '))[2]
      SPECIES=paste(gsub('uncultured ','',SPECIES),'sp.')
    }
    
    if(identical(grep('unclassified',GENUS),integer(0))==F){
      GENUS=unlist(strsplit(SPECIES,split=' '))[2]
      SPECIES=paste(gsub('unclassified ','',SPECIES),'sp.')
    }
    
    SPECIES=gsub(']','',(gsub('[','',SPECIES,fixed = T)),fixed = T)
    GENUS=gsub(']','',(gsub('[','',GENUS,fixed = T)),fixed = T)
    
    if(t==1){NCBI.Tax.backup=NULL}
    
    if(t>1&SPECIES%in%NCBI.Tax.backup$Species){
      tmp=data.frame(
        seq.ID=species$seq.ID[t],
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
      
      search=taxize::classification(SPECIES, db = 'ncbi')
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
        seq.ID=species$seq.ID[t],
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
        
        search=taxize::classification(new.spe.name, db = 'ncbi')
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
          search=taxize::classification(new.spe.name, db = 'ncbi')
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
          seq.ID=species$seq.ID[t],
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
        search=taxize::classification(GENUS, db = 'ncbi')
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
          seq.ID=species$seq.ID[t],
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
      seq.ID=species$seq.ID[t],
      Tax.NBC=paste(
        paste('k__',DO,sep=''),
        paste('p__',PH,sep=''),
        paste('c__',CL,sep=''),
        paste('o__',OR,sep=''),
        paste('f__',FA,sep=''),
        paste('g__',GENUS,sep=''),
        paste('s__',SPECIES,sep=''),sep='; '))
    
    tmp.BLCA=data.frame(
      seq.ID=species$seq.ID[t],
      Tax.BLCA=paste(
        paste('species:',DO,sep=''),
        paste('genus:',PH,sep=''),
        paste('family:',CL,sep=''),
        paste('order:',OR,sep=''),
        paste('class:',FA,sep=''),
        paste('phylum:',GENUS,sep=''),
        paste('superkingdom:',SPECIES,sep=''),sep=';')) 
    
    if(t==1){NCBI.Tax=tmp;NBC=tmp.NBC;BLCA=tmp.BLCA}else{NCBI.Tax=rbind(NCBI.Tax,tmp);NBC=rbind(NBC,tmp.NBC);BLCA=rbind(BLCA,tmp.BLCA)}
    
    NCBI.Tax.backup=NCBI.Tax
  }
  
  write.csv(NCBI.Tax,paste(Target,database,'full_ORF_Taxonomy.csv',sep='_'),row.names = F)
  write.table(NBC,paste(Target,database,'full_ORF_Taxonomy_NBC.txt',sep='_'),col.names = F,row.names = F,quote = F,sep='\t')
  write.table(BLCA,paste(Target,database,'full_ORF_Taxonomy_BLCA.txt',sep='_'),,col.names = F,row.names = F,quote = F,sep='\t')
  print('Taxonomy files generated')
}