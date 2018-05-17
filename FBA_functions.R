###FBA functions

#### Contents
### 1) Contribution/Correlation Analysis functions
### 2) Running MIMOSA on FBA functions
### 3) Plotting functions
### 4) FBA processing functions, misc, network manipulation

############ 1) Contribution/Correlation Analysis Functions ####################

############ Calculate true variance contributions
getContributions = function(met_fluxes_final, spec_codes, kegg_translate, path_key){
  true_met_fluxes = merge(met_fluxes_final, spec_codes, by="Species", all = T)
  if("medium" %in% names(true_met_fluxes)) setnames(true_met_fluxes, "medium", "compound")
  met_fluxes_fill = dcast(true_met_fluxes, compound+SimRun~Species, value.var = "cumulFlux", fun.aggregate=sum)

  var_shares = rbindlist(lapply(spec_codes[,Species], function(y){
    all1 = rbindlist(lapply(spec_codes[,Species], function(x){
      foo = met_fluxes_fill[,cov(get(x), get(y), use="complete.obs"), by=compound]
      foo[,Species:=x]
      return(foo)
    }))
    all1[,Species2:=y]
  }))
  var_shares = var_shares[,sum(V1),by=list(compound, Species)]
  var_shares = merge(var_shares, spec_codes[,list(Code,Species)], by = "Species")
  
  tot_met_fluxes = true_met_fluxes[,sum(cumulFlux), by = list(compound, SimRun)]
  true_met_var = tot_met_fluxes[,list(var(V1), mean(V1)), by = compound]
  setnames(true_met_var, c("V1", "V2"), c("TrueVar", "Mean"))
  var_shares = merge(var_shares, true_met_var, by="compound")
  var_shares[,VarShare:=V1/TrueVar]
  #Merge extra info
  shap_contribs = merge(var_shares, kegg_translate, by.x="compound", by.y = "medium", all.x=T)
  shap_contribs[,niceLab:=gsub("_e0","",Metabolite)]
  shap_contribs = shap_contribs[Metabolite %in% kegg_translate[,Metabolite]]
  met_order = unique(shap_contribs[,list(Metabolite, niceLab, TrueVar)])[order(TrueVar, decreasing=T), niceLab]
  shap_contribs[,niceLab:=factor(niceLab, levels = met_order)]
  # if("KEGG" %in% names(shap_contribs)){
  #   shap_contribs = merge(shap_contribs, path_key, by.x = "KEGG", by.y = "compound", all.x=T, all.y=F)
  #   shap_contribs[is.na(SuperPath), SuperPath:="Other"]
  # }
  shap_contribs = merge(shap_contribs, spec_codes, by="Species", all.x = T)
  return(shap_contribs)
}

############### Calculate correlations and merge with contributions
get_correlation_contrib_comparison = function(fake_spec, fake_mets_melt, shap_contribs, outcome_var = "BinaryContribPos", contrib_threshold = 0.1, spec_codes = make_spec_codes()[1:10]){
  if(!("Sample" %in% names(fake_spec))){
    fake_spec[,Sample:=paste0("run", SimRun, "_", paste0(spec_codes[,Code], collapse=""), "_TP", TimePoint)]
  }
  if(!("Sample" %in% names(fake_mets_melt))){
    fake_mets_melt[,Sample:=paste0("run", SimRun, "_", paste0(spec_codes[,Code], collapse=""), "_TP", TimePoint)]
  }
  if(!("PosVarShare" %in% names(shap_contribs) & "VarShareMagnitude" %in% names(shap_contribs))){
    shap_contribs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0]),0), by=compound]
    shap_contribs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=compound]
  }
  bad_mets = fake_mets_melt[,var(value),by=compound][V1==0|is.na(V1),compound]
  spec_met_corrs = basic_correlation_matrix(fake_spec[,list(Species,Sample,value)], fake_mets_melt[!compound %in% bad_mets,list( compound,Sample, value)], method="spearman")
  
  spec_met_corrs = merge(shap_contribs[!compound %in% bad_mets], spec_met_corrs, by = c("Species", "compound"), all = T)
  spec_met_corrs = spec_met_corrs[!Species %in% c("Inflow", "Outflow")]
  
  spec_met_corrs[,CorTransform:=log((1+estimate)/(1-estimate))/2]
  spec_met_corrs[,CorTransform:=ifelse(estimate==1, max(CorTransform[estimate != 1], na.rm=T)*1.1, CorTransform)]
  spec_met_corrs[,CorTransform:=ifelse(estimate==-1, min(CorTransform[estimate != -1], na.rm=T)*1.1, CorTransform)]
  spec_met_corrs[,CorScale:=(CorTransform-min(CorTransform, na.rm=T))/(max(CorTransform, na.rm=T)-min(CorTransform, na.rm=T))]
  spec_met_corrs[is.na(CorScale), CorScale:=0]
  
  spec_met_corrs[,CorScaleAbs:=(abs(CorTransform)-min(abs(CorTransform), na.rm=T))/(max(abs(CorTransform), na.rm=T)-min(abs(CorTransform), na.rm=T))]
  if(outcome_var == "BinaryContribPos"){
    spec_met_corrs[,BinaryContrib:=ifelse(abs(VarShare) > 0.45,1,0)]
    spec_met_corrs[,BinaryContrib2:=ifelse(VarShare > 0.45,1,0)]
    spec_met_corrs[,BinaryContribPos:=ifelse(PosVarShare > contrib_threshold, 1, 0)]
    spec_met_corrs[,BinaryContribMag:=ifelse(VarShareMagnitude > contrib_threshold*2, 1, 0)]
  }
  spec_met_corrs[,TruePos:=ifelse(p.value < 0.01 & get(outcome_var) == 1,1,0)]
  spec_met_corrs[,FalsePos:=ifelse(p.value < 0.01 & get(outcome_var) == 0,1,0)]
  spec_met_corrs[,TrueNeg:=ifelse(p.value > 0.01 & get(outcome_var) == 0, 1,0)]
  spec_met_corrs[,FalseNeg:=ifelse(p.value > 0.01 & get(outcome_var) == 1, 1,0)]
  spec_met_corrs[,Outcome:=ifelse(TruePos==1,"TruePos", "FalsePos")]
  spec_met_corrs[,Outcome:=ifelse(TrueNeg==1,"TrueNeg", Outcome)]
  spec_met_corrs[,Outcome:=ifelse(FalseNeg==1,"FalseNeg", Outcome)]
  return(spec_met_corrs)
}

########## Get scaled contribution values. This assumes you have all contributors (e.g. also environment)
getScaledContribs = function(spec_met_corrs, contrib_threshold = 0.1){
  if(!("PosVarShare" %in% names(spec_met_corrs) & "VarShareMagnitude" %in% names(spec_met_corrs))){
    spec_met_corrs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0]),0), by=niceLab]
    spec_met_corrs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=niceLab]
    spec_met_corrs[,BinaryContribPos:=ifelse(PosVarShare > contrib_threshold, 1, 0)]
    spec_met_corrs[,BinaryContribMag:=ifelse(VarShareMagnitude > contrib_threshold*2, 1, 0)]
  }
  return(spec_met_corrs)
}


################# Run permutation/subset Shapley contribution analysis
run_shapley_analysis = function(nperm = 15000, true_met_fluxes, spec_codes, save_perms = F, outdir){
  tot_fluxes = true_met_fluxes[,sum(cumulFlux), by=list(medium, SimRun)]
  true_met_var = tot_fluxes[,var(V1), by=medium]
  setnames(true_met_var, "V1", "TrueVar")
  M = length(unique(spec_codes[,Species]))
  
  allCumulMetVars = data.table()
  for(perm_id in 1:nperm){
    cumulMetVars = true_met_var
    spec_order = sample(1:10)
    for(j in 1:length(R1)){
      if(j < length(R1)){
        met_fluxes = true_met_fluxes[!Species %in% spec_codes[spec_order[1:j], Species]]
        #Fill in 0s for all species even if they don't do anything for that compound
        #When they are equal this works fine
        cumulMetFluxes = met_fluxes[,sum(cumulFlux), by = list(medium, SimRun)]
        #calculate cumulative flux var/sd under permutation
        cumulMetVar = cumulMetFluxes[,var(V1), by = medium]
        cumulMetVars = merge(cumulMetVars, cumulMetVar, by = "medium", all.x = T)
        setnames(cumulMetVars, "V1", spec_codes[spec_order[j], Code])
      } else {
        cumulMetVars[,(spec_codes[spec_order[j],Code]):=0]
      }
      if(j > 1){
        cumulMetVars[,paste0("Marg_", spec_codes[spec_order[j], Code]):=get(spec_codes[spec_order[j-1], Code]) - get(spec_codes[spec_order[j], Code])]
      } else {
        cumulMetVars[,paste0("Marg_", spec_codes[spec_order[j], Code]):=TrueVar - get(spec_codes[spec_order[j], Code])]
      }
    }
    cumulMetVars[,OrderID:=perm_id]
    cumulMetVars = cumulMetVars[,c("medium","TrueVar", sort(names(cumulMetVars)[3:ncol(cumulMetVars)])), with=F]
    allCumulMetVars = rbind(allCumulMetVars, cumulMetVars, fill = T)
  }
  if(save_perms) write.table(allCumulMetVars, file = paste0(outdir, "AllCumulMetVars.txt"), quote=F, row.names = F, sep = "\t", col.names = T) #ifelse(perm_id==1 & sub_perm_id==1, T,F))
  allContribs = allCumulMetVars[,lapply(.SD, mean), by=medium, .SDcols = paste0("Marg_", spec_codes[,Code])]
  setnames(allContribs, gsub("Marg_", "", names(allContribs)))
  allContribs = melt(allContribs, variable.name = "Code")
  write.table(allContribs, file = paste0(outdir, "AllShapContribs3.txt"), quote=F, row.names = F, sep = "\t")
}

###################### 2) Running MIMOSA on FBA Functions #################################################

#################### Get reaction content for each sample
FBA_picrust = function(speciesAbunds, specRxnTable, spec_codes){
  all_rxns = specRxnTable[,RxnID]
  all_samps = names(speciesAbunds)[names(speciesAbunds) != "Species"]
  all_spec = sort(speciesAbunds[,Species])
  speciesAbunds = speciesAbunds[,lapply(.SD, as.numeric), .SDcols = all_samps]
  all_codes = spec_codes[order(Species)][Species %in% all_spec, Code]
  specRxnTable = specRxnTable[,lapply(.SD, as.numeric), .SDcols = all_codes]
  rxn_abunds = t(t(as.matrix(speciesAbunds[,all_samps,with=F]))%*%t(as.matrix(specRxnTable[,all_codes,with=F])))
  rxn_abunds = data.table(rxn_abunds, all_rxns)
  setnames(rxn_abunds, c(all_samps, "KO"))
  return(rxn_abunds)
}

#################### Get reaction content for each species
FBA_picrust_contribs = function(speciesAbunds, specRxnTable_long){
  all_rxns = specRxnTable_long[,unique(RxnID)]
  all_samps = names(speciesAbunds)[names(speciesAbunds) != "Species"]
  all_spec = speciesAbunds[,Species]
  speciesAbunds = speciesAbunds[,lapply(.SD, as.numeric), .SDcols = all_samps]
  speciesAbunds[,Species:=all_spec]
  spec_abunds = melt(speciesAbunds, id.var = "Species", variable.name = "Sample")
  all_codes = specRxnTable_long[,unique(Code)]
  #specRxnTable = specRxnTable[,lapply(.SD, as.numeric), .SDcols = all_codes]
  contribs = merge(specRxnTable_long, spec_abunds, by = "Species", allow.cartesian = T)
  contribs[,singleMusicc:=CopyNum*value]
  return(contribs)
}


########## Take list of FBA models and get components needed for MIMOSA
# remove_rev: whether to remove reversible reactions from S mat
# missing_rxns: whether to include reactions with no annotated genes (gap-filling, transport rxns)
build_model_components = function(all_mods, remove_rev = T, missing_rxns = F, spec_codes, agora = F){
  #ID reversible reactions
  reversible_rxns = data.table()
  for(j in 1:length(all_mods)){
    reversible_rxns = rbind(reversible_rxns, data.table(Rxn = unlist(all_mods[[j]]$rxns), Rev = unlist(all_mods[[j]]$rev), Code = spec_codes[j,Code], Species = spec_codes[j, Species]))
  }
  reversible_rxns = reversible_rxns[Rxn != "biomass0"]
  reversible_rxns = reversible_rxns[!grepl("EX_", Rxn)]
  
  #get S matrices
  all_S_mat = list()
  for(j in 1:length(all_mods)){
    all_S_mat[[j]] = as.matrix(all_mods[[j]]$S)
    row.names(all_S_mat[[j]]) = unlist(all_mods[[j]]$mets)
    colnames(all_S_mat[[j]]) = unlist(all_mods[[j]]$rxns)
  }
  
  all_comps = lapply(all_S_mat, function(x){ return(row.names(x))})
  #   length(Reduce(intersect, all_comps)) #572
  #   length(Reduce(union, all_comps)) #1448
  
  all_rxns = lapply(all_S_mat, function(x){ return(colnames(x))})
  #   length(Reduce(intersect, all_rxns)) #389
  #   length(Reduce(union, all_rxns)) #1930
  #Well, lots of differences
  
  all_S_mat = lapply(all_S_mat, function(x){ 
    foo = data.table(Compound = row.names(x), x)
    foo[,biomass0:=NULL] #don't need this
    foo = foo[,which(!grepl("EX_",names(foo))),with=F]
    return(foo)})
  
  all_S_mats = all_S_mat[[1]]
  for(j in 2:length(all_S_mat)){
    all_S_mats = merge(all_S_mats, all_S_mat[[j]], by=intersect(names(all_S_mats), names(all_S_mat[[j]])), all=T)
  }
  #Community EMM
  #Get rid of duplicates
  all_S_mats = all_S_mats[,lapply(.SD, function(x){ if(length(x[!is.na(x)]) > 0 ) return(unique(x[!is.na(x)])) else return(0) }), by=Compound]
  setkey(all_S_mats, NULL)
  
  if(agora == F){
    #Get counts of each Rxn for every species
    all_geneRxn_mats = lapply(1:length(all_mods), function(x){ 
      mat = data.table(t(as.matrix(all_mods[[x]]$rxnGeneMat)))
      setnames(mat, unlist(all_mods[[x]]$rxns))
      mat[,Gene:=unlist(all_mods[[x]]$genes)]
      mat[,Species:=spec_codes[x,Code]]
      mat = mat[Gene != "kb"]
      mat[Gene=="Unknown", Gene:=paste0("Unknown",spec_codes[x,Code])]
      return(mat)
    })
    if(missing_rxns){
      #For every reaction annotated in a species' S mat but without any gene in the gene mat, add unknown gene
      #Mostly transport reactions
      for(k in 1:length(all_mods)){
        missing_genes = names(all_S_mat[[k]])[which(unlist(all_geneRxn_mats[[k]][,lapply(.SD, function(x){ sum(x) }), .SDcols = names(all_S_mat[[k]])[which(names(all_S_mat[[k]]) != "Compound")]])==0)+1]
        for(i in 1:length(missing_genes)){
          set(all_geneRxn_mats[[k]], all_geneRxn_mats[[k]][,grep("Unknown",Gene)[1]], which(names(all_geneRxn_mats[[k]])==missing_genes[i]), 1)
        }
      }
    }
    all_geneRxn_mat = rbindlist(all_geneRxn_mats, fill = T)
    for (j in 1:ncol(all_geneRxn_mat)){
      set(all_geneRxn_mat,which(is.na(all_geneRxn_mat[[j]])),j,0)
    }
    specRxns = all_geneRxn_mat[,lapply(.SD, sum), by=Species, .SDcols=which(!names(all_geneRxn_mat)  %in% c("Gene","Species", "biomass0") & !grepl("EX_", names(all_geneRxn_mat)))] #
    specRxns = specRxns[,c(1, which(colSums(specRxns[,2:ncol(specRxns),with=F])!=0)+1),with=F]
    specRxns = data.table(t(specRxns[,2:ncol(specRxns),with=F]), RxnID = names(specRxns)[names(specRxns) != "Species"])
    setnames(specRxns, c(spec_codes[,Code], "RxnID"))
  } else {
    specRxns = rbindlist(lapply(1:length(all_mods), function(x){
      rxn_dat = data.table(Code = spec_codes[x,Code], RxnID = unlist(all_mods[[x]]$rxns), value = 1)
    }))
    specRxns = dcast(specRxns, RxnID~Code, value.var="value", fill = 0)
  }
  
  #Separate reactions annotated differently in different species
  dups = unique(all_S_mats[duplicated(Compound),Compound])
  if(length(dups) > 0){
    conflict_rxns = list()
    for(i in 1:length(dups)){
      nonzero = all_S_mats[Compound==dups[i]][,names(which(sapply(.SD, uniqueN)!=1))]
      conflict_rxns[[i]] = names(all_S_mats[Compound==dups[i]][,nonzero,with=F])
    }
    all_conflict_rxns = sort(unique(unlist(conflict_rxns)))
    
    for(j in 1:length(all_conflict_rxns)){
      dup_vals = all_S_mats[Compound %in% dups, unique(get(all_conflict_rxns[j])),by=Compound][duplicated(Compound),Compound]
      vals = all_S_mats[Compound %in% dup_vals, list(Compound,get(all_conflict_rxns[j]))]
      #We need to figure out which compound coefficients go with which from the original S mat
      all_options = list()
      spec_list = list()
      count_opt = 1
      any_match = F
      for(m in 1:length(spec_codes[,Species])){ #check which version each species has
        if(all_conflict_rxns[j] %in% names(all_S_mat[[m]])){ 
          opt1 = all_S_mat[[m]][Compound %in% dup_vals,list(Compound,get(all_conflict_rxns[j]))][order(Compound)]
          any_match = F
          q = 1
          while(any_match ==F & q <= length(all_options)){
            if(all(opt1 == all_options[[q]])){
              spec_list[[q]] = c(spec_list[[q]], spec_codes[m,Species])
              any_match = T
            } 
            q = q+1
          }
          if(any_match == F){ #if no matching version for this yet
            all_options[[count_opt]] = opt1
            spec_list[[count_opt]] = spec_codes[m,Species]
            count_opt = count_opt + 1
          }
        }
      }
      for(k in 1:length(all_options)){
        setnames(all_options[[k]], "V2", paste0(all_conflict_rxns[j], "_",k))
        all_S_mats = merge(all_S_mats, all_options[[k]], by="Compound", all = T)
        set(all_S_mats, i=which(is.na(all_S_mats[,get(paste0(all_conflict_rxns[j], "_",k))])), j=which(names(all_S_mats)==paste0(all_conflict_rxns[j], "_",k)), value = all_S_mats[is.na(get(paste0(all_conflict_rxns[j], "_",k))), get(all_conflict_rxns[j])]) 
        specRxns = rbind(specRxns, specRxns[RxnID==all_conflict_rxns[j]])
        specRxns[nrow(specRxns), RxnID:=paste0(RxnID, "_", k)]
        spec_codes_bad = spec_codes[(Species %in% unlist(spec_list)) & !(Species %in% spec_list[[k]]),Code]
        specRxns[nrow(specRxns), (spec_codes_bad):=0] 
        reversible_rxns[Rxn==all_conflict_rxns[j] & Species %in% spec_list[[k]], Rxn:=paste0(all_conflict_rxns[j], "_",k)]
        #rename rxns in the original S_mats
        spec_ids = match(spec_list[[k]], spec_codes[,Species])
        for(m in 1:length(spec_ids)){
          setnames(all_S_mat[[spec_ids[m]]], all_conflict_rxns[j], paste0(all_conflict_rxns[j], "_", k))
        }
      }
      all_S_mats[,(all_conflict_rxns[j]):=NULL]
      specRxns = specRxns[RxnID != all_conflict_rxns[j]]
      setkey(all_S_mats, NULL)
    }
    all_S_mats = unique(all_S_mats)
  }
  
  #Merge specRxns with reversibility info
  specRxns_melt = melt(specRxns, id.var = "RxnID", variable.name = "Code", value.name = "CopyNum")
  specRxns_melt = merge(specRxns_melt, reversible_rxns, by.x = c("RxnID","Code"), by.y = c("Rxn", "Code"), all.x=T)
  specRxns_melt = specRxns_melt[CopyNum != 0]
  
  if(remove_rev){
    specRxns_melt_rev = copy(specRxns_melt)
    specRxns_rev = copy(specRxns) ##So that we can modify one and not the other later
    specRxns_melt = specRxns_melt[Rev.V1==0]
    rev_rxn_ids = specRxns_melt_rev[Rev.V1==1, unique(RxnID)] #Reactions which are reversible for any subset of species
    for(k in 1:length(rev_rxn_ids)){ #For species that have the reversible version, remove from specRxns
      spec_fix = specRxns_melt_rev[RxnID==rev_rxn_ids[k] & Rev.V1==1, Code]
      set(specRxns, i = specRxns[,which(RxnID==rev_rxn_ids[k])], j=which(names(specRxns) %in% spec_fix), 0)
    }
    specRxns = specRxns[rowSums(specRxns[,which(names(specRxns) != "RxnID"),with=F]) != 0]
    
    all_S_mats_rev = all_S_mats
    all_S_mats = all_S_mats[,names(all_S_mats) %in% c("Compound", specRxns[,unique(RxnID)]),with=F]
    all_S_mat_rev = all_S_mat
    for(j in 1:length(all_S_mat)){
      all_S_mat[[j]] = all_S_mat[[j]][,names(all_S_mat[[j]]) %in% c("Compound", specRxns_melt[Species==spec_codes[j,Species],RxnID]),with=F]
    }
  }
  return(list(all_S_mat, all_S_mats, specRxns, specRxns_melt, reversible_rxns, all_S_mat_rev, all_S_mats_rev, specRxns_rev, specRxns_melt_rev))
}

############ Make reaction IDs in flux data consistent with mimosa_pieces
fix_rxn_ids = function(run_data, specRxns_melt_rev){
  changed_rxns = run_data$fluxes[!Rxn %in% specRxns_melt_rev[,unique(RxnID)], unique(Rxn)]
  if(length(changed_rxns) > 0){
    for(i in 1:length(changed_rxns)){
      replacements = specRxns_melt_rev[grepl(changed_rxns[i], RxnID)] #core needs to be the same
      replacements[,Rxn:=gsub("_[1-2]$","", RxnID)]
      spec_list = replacements[,unique(Species)]
      if(length(replacements[,unique(Rxn)]) > 1) stop("Something weird, multiple core reactions")
      if(nrow(replacements) > 0){
        if(replacements[,unique(Rxn)]==changed_rxns[i]){
          run_data$fluxes = merge(run_data$fluxes, replacements[,list(RxnID, Species, Rxn)], by=c("Rxn", "Species"), all.x = T, all.y = F)
          run_data$fluxes[!is.na(RxnID),Rxn:=RxnID]
          run_data$fluxes[,RxnID:=NULL]
          run_data$glpk_fluxes = merge(run_data$glpk_fluxes, replacements[,list(RxnID, Species, Rxn)], by=c("Rxn", "Species"), all.x = T, all.y = F)
          run_data$glpk_fluxes[!is.na(RxnID),Rxn:=RxnID]
          run_data$glpk_fluxes[,RxnID:=NULL]
        }
      }##Don't do anything if replacement does not match or is not there
    }
  }
  return(run_data)
}

################ Set reversible reactions in the community model as irreversible if they go the same way at 90% of time points
fix_rev_rxns = function(ten_spec_vary, old_mimosa_pieces){
  mimosa_pieces = copy(old_mimosa_pieces) #just in case, don't want to directly modify these
  #compile all fluxes
  all_fluxes = rbindlist(lapply(1:length(ten_spec_vary), function(x){
    flux_sub = ten_spec_vary[[x]]$fluxes
    flux_sub[,SimRun:=x]
    return(flux_sub)
  }))
  rev_rxns = mimosa_pieces[[5]]
  #Get rev rxns
  rev_rxn_ids = rev_rxns[Rev.V1==1]
  #Look at only fluxes of rev rxns
  all_fluxes = all_fluxes[!grepl("EX",Rxn)] #going to remove transport reactions anyway
  total_tps = all_fluxes[,length(value),by=list(Species,Rxn)][1,V1]
  confint_upper = total_tps*0.95
  #Remove rxns that only ever have zero flux
  all_fluxes[,BadRxn:=any(value != 0), by=list(Species,Rxn)]
  all_fluxes = all_fluxes[BadRxn==T]
  all_fluxes = merge(all_fluxes, rev_rxn_ids, by = c("Rxn","Species"), all.x = F, all.y = T)
  
  #Get reversible reactions to make irreversible
  rxn_counts = all_fluxes[,list(length(cumulFlux[cumulFlux <= 0]), length(cumulFlux[cumulFlux >= 0])), by=list(Species,Rxn)]
  #A reaction can be above confint_upper for both if it has a lot of zeros - remove these
  rxn_counts = rxn_counts[!(V2 > confint_upper & V1 > confint_upper)]
  new_pos_rxns = rxn_counts[V2 > confint_upper]
  new_pos_rxns[,RxnChange:=1]
  new_neg_rxns = rxn_counts[V1 > confint_upper]
  new_neg_rxns[,RxnChange:=-1]
  new_rxns = rbind(new_pos_rxns, new_neg_rxns)
  new_rxns[,RxnTypeCount:=length(unique((RxnChange))),by=Rxn]
  #Split up rxns that go different ways in different species
  switch_rxns = new_rxns[RxnTypeCount==2,unique(Rxn)]
  new_rxns[,NewRxn:=Rxn]
  new_rxns[Rxn %in% switch_rxns & RxnChange==1, NewRxn:=paste0(Rxn,"_1")]
  new_rxns[Rxn %in% switch_rxns & RxnChange==-1, NewRxn:=paste0(Rxn,"_2")]
  
  #Now fix in every mimosa piece  
  mimosa_pieces[[5]] = merge(mimosa_pieces[[5]], new_rxns[,list(Species,Rxn,NewRxn, RxnChange)], all.x=T, all.y = F, by=c("Species","Rxn"))
  mimosa_pieces[[5]][!is.na(RxnChange), Rev.V1:=0]
  mimosa_pieces[[5]][!is.na(RxnChange), Rxn:=NewRxn]
  
  #replace old IDs with split up IDs
  new_ten_spec_vary = list()
  for(j in 1:length(ten_spec_vary)){
    #remove transporters
    y = copy(ten_spec_vary[[j]])
    y$fluxes = y$fluxes[!grepl("EX", Rxn)]
    y$fluxes = merge(y$fluxes, new_rxns[,list(Rxn, Species, NewRxn, RxnChange)], by = c("Rxn", "Species"), all.x = T, all.y = F)
    y$fluxes[!is.na(NewRxn) & Rxn != NewRxn, Rxn:=NewRxn]
    y$fluxes[,NewRxn:=NULL]
    #Make fluxes positive since we will switch reaction to opposite direction
    y$fluxes[RxnChange==-1, value:=-1*value]
    y$fluxes[RxnChange==-1, totalFlux:=-1*totalFlux]
    y$fluxes[RxnChange==-1, cumulFlux:=-1*cumulFlux]
    y$fluxes[,RxnChange:=NULL]
    new_ten_spec_vary = append(new_ten_spec_vary, list(y))
  }
  ten_spec_vary = new_ten_spec_vary
  
  spec_codes[,SpecNum:=1:10]
  new_rxns = merge(new_rxns, spec_codes, by="Species", all.x=T)
  
  #Fix other mimosa pieces
  melted_pieces1 = lapply(mimosa_pieces[[1]], melt, variable.name = "Rxn")
  melted_pieces6 = lapply(mimosa_pieces[[6]], melt, variable.name = "Rxn")
  
  #Need to add in new non-reversible rxns for part 1, they will already be there for part 6
  melted_pieces1 = lapply(1:length(melted_pieces1), function(x){
    spec = spec_codes[x,Species]
    new_pieces = melted_pieces6[[x]][Rxn %in% new_rxns[Species==spec,Rxn] & !Rxn %in% melted_pieces1[[x]][,Rxn]] #new non-rev rxns
    foo = rbind(melted_pieces1[[x]], new_pieces)
    foo = merge(foo, new_rxns[Species == spec,list(Rxn,NewRxn, RxnChange)], by = "Rxn", all.x = T, all.y = F)
    foo[Rxn != NewRxn, Rxn:=NewRxn]
    foo[RxnChange == -1, value:=-1*value]
    foo[,NewRxn:=NULL]
    foo[,RxnChange:=NULL]
    return(foo)
  })
  
  melted_pieces6 = lapply(1:length(melted_pieces6), function(x){
    spec = spec_codes[x,Species]
    foo = merge(melted_pieces6[[x]], new_rxns[Species == spec,list(Rxn,NewRxn, RxnChange)], by = "Rxn", all.x = T, all.y = F)
    foo[Rxn != NewRxn, Rxn:=NewRxn]
    foo[RxnChange == -1, value:=-1*value]
    foo[,NewRxn:=NULL]
    foo[,RxnChange:=NULL]
    return(foo)
  })
  
  for(j in 1:nrow(new_rxns)){
    #if old reaction is still there, save info
    if(new_rxns[j,Rxn] %in% names(mimosa_pieces[[7]])){
      old_rxn_info = mimosa_pieces[[7]][,get(new_rxns[j,Rxn])]
      if(new_rxns[j,NewRxn != Rxn]){
        mimosa_pieces[[7]][,(new_rxns[j,Rxn]) :=NULL]
      }
    } else if(new_rxns[j,NewRxn] %in% names(mimosa_pieces[[7]])){
      old_rxn_info = mimosa_pieces[[7]][,get(new_rxns[j,NewRxn])] #Already fixed it
    } else{ #Changed rxn ID
      matching_rxn = names(mimosa_pieces[[7]])[grepl(new_rxns[j,Rxn], names(mimosa_pieces[[7]]), fixed = T)]
      old_rxn_info = mimosa_pieces[[7]][,get(matching_rxn)]
    }
    
    if(new_rxns[j,RxnChange==1]){ #If forward
      
      if(!(new_rxns[j,NewRxn] %in% names(mimosa_pieces[[2]]))){ #If we haven't gotten to this with another species
        mimosa_pieces[[2]][,new_rxns[j,NewRxn]:=old_rxn_info] 
        mimosa_pieces[[7]][,new_rxns[j,NewRxn]:=old_rxn_info] 
      }
    } else { #Reverse
      if(!(new_rxns[j,NewRxn] %in% names(mimosa_pieces[[2]]))){ #If we haven't gotten to this with another species
        mimosa_pieces[[2]][,new_rxns[j,NewRxn]:=-1*old_rxn_info] 
        mimosa_pieces[[7]][,new_rxns[j,NewRxn]:=-1*old_rxn_info]
      }
    }
    }
    for(k in 1:length(melted_pieces1)){
      mimosa_pieces[[1]][[k]] = dcast(melted_pieces1[[k]], Compound~Rxn, value.var = "value", fun.aggregate = sum)
      mimosa_pieces[[6]][[k]] = dcast(melted_pieces6[[k]], Compound~Rxn, value.var = "value", fun.aggregate = sum)
    }
    mimosa_pieces[[6]] = lapply(mimosa_pieces[[6]], function(x){
      return(x[,which(!grepl("EX_",names(x))), with = F])
    })
    #Remove transport reactions
    mimosa_pieces[[7]] = mimosa_pieces[[7]][,which(!grepl("EX_", names(mimosa_pieces[[7]]))), with = F]
    mimosa_pieces[[8]] = mimosa_pieces[[8]][!grepl("EX_", RxnID)]
    mimosa_pieces[[9]] = mimosa_pieces[[9]][!grepl("EX_", RxnID)]
    #Add additional nonrev rxns
    new_specRxns = merge(new_rxns[,list(Rxn,Species,RxnChange, NewRxn)], mimosa_pieces[[9]], all.y = F, by.x = c("Rxn", "Species"), by.y = c("RxnID", "Species"))
    new_specRxns[,Rev.V1:=0]
    setnames(new_specRxns, "NewRxn", "RxnID")
    
    mimosa_pieces[[4]] = rbind(mimosa_pieces[[4]], new_specRxns[,list(RxnID, Code,CopyNum, Rev.V1, Species)])
    #Fix remaining mimosa pieces
    mimosa_pieces[[9]] = merge(mimosa_pieces[[9]], new_rxns, by.x = c("RxnID", "Species", "Code"), by.y = c("Rxn", "Species", "Code"), all = T)
    mimosa_pieces[[9]][!is.na(RxnChange), Rev.V1:=0] #set for merged rxns
    mimosa_pieces[[9]][RxnID != NewRxn, RxnID:=NewRxn]
    mimosa_pieces[[9]] = mimosa_pieces[[9]][,list(RxnID, Code, CopyNum, Rev.V1, Species)]
    #Remove renamed/replaced rxns
    rxns_to_remove = new_rxns[Rxn != NewRxn, Rxn]
    mimosa_pieces[[9]] = mimosa_pieces[[9]][!RxnID %in% rxns_to_remove]
    
    mimosa_pieces[[8]] = dcast(mimosa_pieces[[9]], RxnID~Code, value.var = "CopyNum")
    mimosa_pieces[[8]][is.na(mimosa_pieces[[8]])] = 0
    mimosa_pieces[[8]] = mimosa_pieces[[8]][,c(spec_codes[,Code], "RxnID"), with=F]
    mimosa_pieces[[3]] = dcast(mimosa_pieces[[4]], RxnID~Code, value.var = "CopyNum")
    mimosa_pieces[[3]][is.na(mimosa_pieces[[3]])] = 0
    mimosa_pieces[[3]] = mimosa_pieces[[3]][,c(spec_codes[,Code], "RxnID"), with=F]
    return(list(ten_spec_vary, mimosa_pieces))
}
  
  #################### Run core MIMOSA analysis on simulation dataset
  #gene_type either species, totalFlux, cumulFlux, or genes
run_all_metabolites_FBA = function(run_prefix, fake_spec = "", fake_genes = "", fake_mets, gene_type = "species", species_rxns, rxn_mets, mimosa_pieces = "", spec_codes, degree_filter = 40, cor_method = "spearman", correction = "fdr", nperm = 10000, nonzero_filter = 4, all_fluxes = "", kegg_translate = ""){
    if(gene_type == "species"){
      if(identical(fake_spec, "")) stop("Include species to use for gene abundances")
      #Using only nonrev rxns for consistency/safety
      if(mimosa_pieces != ""){
        specRxns = mimosa_pieces[[3]][,lapply(.SD, as.numeric), .SDcols = spec_codes[,Code]]
        specRxns[,RxnID:=mimosa_pieces[[3]][,RxnID]]
      } else {
        specRxns = species_rxns
      }
      fake_rxn_abunds = FBA_picrust(fake_spec, specRxns, spec_codes)
    } else if(gene_type == "totalFlux" | gene_type == "cumulFlux"){
      if(identical(all_fluxes, "")) stop("Include all_fluxes to use reaction fluxes")
      #all_fluxes[,Sample:=paste0(run_prefix, SimRun, "_", all_spec, "_TP", TimePoint)]
      fake_rxn_abunds = dcast(all_fluxes, Rxn~Sample, value.var = "cumulFlux", fun.aggregate = sum) #add over all species
      setnames(fake_rxn_abunds, "Rxn", "KO")
    } else {
      if(identical(fake_genes, "")) stop("Must include already set up gene abundance matrix")
      fake_rxn_abunds = fake_genes
    }
    
    if(mimosa_pieces != ""){
      #Exclude rev reactions from S mat
      all_S_mats = mimosa_pieces[[2]]
    } else {
      all_S_mats = emm
    }
    all_comps = all_S_mats[,Compound]
    subjects = intersect(names(fake_rxn_abunds), names(fake_mets))
    if(length(subjects) < 2) stop("Sample names not consistent between genes and metabolites")
    
    emm = data.frame(all_S_mats[,2:ncol(all_S_mats),with=F])
    row.names(emm) = gsub("[e]", "[env]", all_comps, fixed=T)
    norm_kos = fake_rxn_abunds
    emm = emm[,names(emm) %in% norm_kos[,KO]]
    cmp_mat = get_cmp_scores(emm, norm_kos)
    
    #get mets
    if(!"medium" %in% names(fake_mets) & !"compound" %in% names(fake_mets)){
      fake_mets = merge(fake_mets, kegg_translate, by="Metabolite", all.x=T, all.y=F)
    }
    if("compound" %in% names(fake_mets)){
      metIDs = fake_mets[,compound] 
      met_name = "compound"
    } else {
      metIDs = fake_mets[,medium]
      met_name = "medium"
    }
    shared_mets = metIDs[metIDs %in% row.names(emm)] 
    setkey(cmp_mat, compound)
    all_comparisons = vector("list",length(shared_mets))
    
    for(j in 1:length(shared_mets)){
      good_subs = intersect(names(fake_mets)[which(!is.na(unlist(fake_mets[shared_mets[j]])))], names(cmp_mat)[which(!is.na(unlist(cmp_mat[shared_mets[j],subjects,with=F])))])
      prmt_vector = unlist(cmp_mat[shared_mets[j],good_subs,with=F])
      met_vector = unlist(fake_mets[shared_mets[j],good_subs,with=F])
      #check for too many 0s
      if(length(met_vector[met_vector!=0 & !is.na(met_vector)]) <= nonzero_filter | length(prmt_vector[prmt_vector!=0]) <= nonzero_filter | length(unique(met_vector)) < 2 | length(unique(prmt_vector)) < 2){
        all_comparisons[[j]] = NA
      }else{
        met_mat = make_pairwise_met_matrix(shared_mets[j], cmp_mat[,c(good_subs, "compound"),with=F])
        metabol_mat = make_pairwise_met_matrix(shared_mets[j], fake_mets[,c(good_subs,met_name),with=F])
        test = mantel_2sided(met_mat,metabol_mat,method=cor_method,permutations = nperm, direction = "pos")
        test_n = mantel_2sided(met_mat,metabol_mat,method=cor_method,permutations = nperm, direction = "neg")
        all_comparisons[[j]] = list(ID = shared_mets[j], PRMT = met_mat, Mets = metabol_mat, Mantel = list(test,test_n))
      }
    }
    failed_mets = shared_mets[which(is.na(all_comparisons))]
    shared_mets = shared_mets[which(!is.na(all_comparisons))]
    all_comparisons = all_comparisons[which(!is.na(all_comparisons))]
    
    correction = "fdr"
    cors_s = sapply(all_comparisons,function(x){return(x$Mantel[[1]]$statistic)})
    pvals_s = sapply(all_comparisons,function(x){return(x$Mantel[[1]]$signif)})
    pvals2_s = correct(pvals_s, method = correction)
    
    cors_n = sapply(all_comparisons,function(x){return(x$Mantel[[2]]$statistic)})
    pvals_n = sapply(all_comparisons,function(x){return(x$Mantel[[2]]$signif)})
    pvals2_n = correct(pvals_n, method = correction)
    
    node_data = data.table::data.table(compound = shared_mets, Correlation = cors_s, PValPos = pvals_s, QValPos = pvals2_s, PValNeg = pvals_n, QValNeg = pvals2_n)
    setkey(node_data,compound)
    
    node_data[,PredictionType:=ifelse(QValPos < 0.1, "Consistent", "Inconsistent")]
    node_data[,PredictionType:=ifelse(QValNeg < 0.1, "Contrasting", PredictionType)]
    
    #save melted version for easier post processing
    cmp_mat_save = melt(cmp_mat, id.var = "compound", variable.name = "Sample", variable.factor = F)
    return(list(all_comparisons, node_data, failed_mets, cmp_mat_save))
}
  
  ############### Get single-species CMPs for simulation data
singleSpecCMP_FBA = function(all_otus, all_koAbunds_byOTU, emms) {
    if (length(all_otus) != length(all_koAbunds_byOTU) | length(all_otus) != length(emms)) 
      stop("Problem! OTU lists don't match")
    cmps_alone = vector("list", length(all_otus))
    for (k in 1:length(all_otus)) {
      all_comps = emms[[k]][,Compound]
      emm = data.frame(emms[[k]][,2:ncol(emms[[k]]),with=F])
      row.names(emm) = all_comps
      #emm = emm[,names(emm) %in% all_koAbunds_byOTU[[k]][,KO]]
      if(!is.null(all_koAbunds_byOTU[[k]])){
        cmps_alone[[k]] = get_cmp_scores(emm, all_koAbunds_byOTU[[k]])
      } else {
        cmps_alone[[k]] = NULL
      }
    }
    return(cmps_alone)
}
  
### Format contributions table, add true contribution information
fix_mimosa_contribs_table = function(mimosa_spec_contribs, node_data_all, shap_contribs_all, id_vars = c("compound", "Species", "Version", "ContribMethod"), comp_list = c(), contrib_threshold_pos = 0.1, contrib_threshold_mag = 0.2){
  #id_vars are ids in mimosa_contribs table
  mimosa_spec_contribs = merge(mimosa_spec_contribs, node_data_all, by = intersect(id_vars, names(node_data_all)),all.x = T)
  mimosa_spec_contribs = merge(mimosa_spec_contribs, shap_contribs_all, by = intersect(id_vars, names(shap_contribs_all)), all.x = T)
  mimosa_spec_contribs[is.na(Pass), Pass:=0]
  if(length(comp_list) > 0) mimosa_spec_contribs = mimosa_spec_contribs[compound %in% comp_list]
  mimosa_spec_contribs[,BinaryContribPos:=ifelse(PosVarShare > contrib_threshold_pos, 1, 0)]
  mimosa_spec_contribs[,BinaryContribMag:=ifelse(VarShareMagnitude > contrib_threshold_mag, 1, 0)]
  mimosa_spec_contribs[,MOutcome:=ifelse(Pass==1 & BinaryContribPos==1, "True positive", "True negative")]
  mimosa_spec_contribs[,MOutcome:=ifelse(Pass==1 & BinaryContribPos==0, "False positive", MOutcome)]
  mimosa_spec_contribs[,MOutcome:=ifelse(Pass==0 & BinaryContribPos==1, "False negative", MOutcome)]
  mimosa_spec_contribs[,PassWords:=ifelse(Pass==1, "MIMOSA\ncontributor", "MIMOSA\nnon-contributor")]
  
  mimosa_spec_contribs[,Cor2:=ifelse(Cor==1, 0.9999999999,Cor)]
  mimosa_spec_contribs[,Cor2:=ifelse(PredictionType != "Consistent", 0, Cor2)] #Only clearly defined for consistent metabolites
  mimosa_spec_contribs[,CorTransform:=log((1+Cor2)/(1-Cor2))/2]
  mimosa_spec_contribs[,CorScale:=(CorTransform-min(CorTransform, na.rm=T))/(max(CorTransform, na.rm=T)-min(CorTransform, na.rm=T)), by=c(id_vars[!id_vars %in% c("Species", "compound")])]
  mimosa_spec_contribs[is.na(CorScale), CorScale:=0]
  return(mimosa_spec_contribs)
}  

mimosa_contrib_comparison = function(mimosa_spec_contribs, shap_contribs, contrib_var = "BinaryContrib"){
  #shap_contribs needs to have BinaryContrib variable already
  mimosa_spec_contribs = merge(mimosa_spec_contribs, shap_contribs, by = c("compound", "Species"))
  mimosa_spec_contribs = mimosa_spec_contribs[compound %in% shap_contribs[,unique(compound)]]
  mimosa_spec_contribs[,MOutcome:=ifelse(Pass==1 & get(contrib_var)==1, "True positive", "True negative")]
  mimosa_spec_contribs[,MOutcome:=ifelse(Pass==1 & get(contrib_var)==0, "False positive", MOutcome)]
  mimosa_spec_contribs[,MOutcome:=ifelse(Pass==0 & get(contrib_var)==1, "False negative", MOutcome)]
  mimosa_spec_contribs[,PassWords:=ifelse(Pass==1, "MIMOSA\ncontributor", "MIMOSA\nnon-contributor")]
  return(mimosa_spec_contribs)
}

################# Get compatible KO data for each species  
get_kegg_genomes_species = function(spec_names, pan_genome = T){ #if pan_genome, will take union of all genomes, otherwise will just return genome of first reference genome
  final_genomes = data.table()
  for(j in 1:length(spec_names)){
    genome_ids = gsub("genome:","", names(keggFind(spec_names[j], database = "genome")))
    if(length(genome_ids)==0){ 
      #stop(paste0("No genomes found for ", spec_names[j])) 
      cat(paste0("No genomes found for ", spec_names[j], "\n"))
      next
    }
    if(pan_genome){
      ko_list = list()
      for(k in 1:length(genome_ids)){
        ko_list[[k]] = data.table(table(gsub("ko:","", keggLink(target = "ko", genome_ids[k])))) #Get copy number info
        ko_list[[k]][,Genome:=genome_ids[k]]
      }
      ko_list = rbindlist(ko_list)
      setnames(ko_list, "V1", "KO")
      ko_list = ko_list[,max(N), by=KO]
      setnames(ko_list, "V1", "CopyNum")
    } else{
      ko_list = data.table(table((gsub("ko:","", keggLink(target = "ko", genome_ids[1])))))
      setnames(ko_list, c("KO", "CopyNum"))
    }
    ko_list[,SpeciesName:=spec_names[j]]
    final_genomes = rbind(final_genomes, ko_list)
  }
  return(final_genomes)
}

################## Process downloaded IMG KO data for each species, to use for MIMOSA with KEGG
process_downloaded_img_data = function(spec_codes, genome_ids_file, gene_annot_dir, pan_genome = T){
  all_genome_info = fread(genome_ids_file, colClasses = c(rep("character", 6), rep("numeric",2)))
  all_genome_info[,SpeciesName:=sapply(`Genome Name / Sample Name`, function(x){ 
    foo = strsplit(x, " ")[[1]][1:2]
    foo[2] = tolower(foo[2])
    return(paste0(foo, collapse = " "))})]
  all_genome_info = all_genome_info[SpeciesName %in% spec_codes[,SpeciesName]] #Removes genomes with non-corresponding species names
  all_genome_info[,DescriptionMatch:=sapply(`Genome Name / Sample Name`, function(x){ any(sapply(spec_codes[,Description], grepl, x = x, ignore.case = T))})]
  missing_species = spec_codes[!SpeciesName %in% all_genome_info[DescriptionMatch==T, SpeciesName], SpeciesName]
  all_genome_info = all_genome_info[DescriptionMatch==T | SpeciesName %in% missing_species]
  all_gene_annots = data.table()
  if(pan_genome){
    for(j in 1:nrow(all_genome_info)){
      gene_annot1 = fread(paste0(gene_annot_dir, all_genome_info[j,`IMG Genome ID`], ".annot.xls"))
      gene_annot1 = gene_annot1[grepl("KO:",Source)]
      gene_annot1[,GenomeID:=all_genome_info[j,`IMG Genome ID`]]
      gene_annot1[,SpeciesName:=all_genome_info[j,SpeciesName]]
      gene_annot1[,Source:=gsub("KO:","", Source)]
      gene_annot1 = gene_annot1[,length(unique(`Locus Tag`)), by=list(Source, GenomeID, SpeciesName)]
      setnames(gene_annot1, c("V1","Source"), c("CopyNum","KO"))
      all_gene_annots = rbind(all_gene_annots, gene_annot1)
    }
    all_gene_annots[,CopyNum:=max(CopyNum), by = list(SpeciesName, KO)]
  } else {
    #Just take the biggest genome?
    genomes_to_use = all_genome_info[,`IMG Genome ID`[which.max(`Gene Count`)], by = SpeciesName]
    for(j in 1:length(genomes_to_use)){
      gene_annot1 = fread(paste0(gene_annot_dir, genomes_to_use[j], ".annot.xls"))
      gene_annot1 = gene_annot1[grepl("KO:",Source)]
      gene_annot1[,GenomeID:=genomes_to_use[j]]
      gene_annot1[,SpeciesName:=all_genome_info[`IMG Genome ID`==genomes_to_use[j],SpeciesName]]
      gene_annot1[,Source:=gsub("KO:","", Source)]
      gene_annot1 = gene_annot1[,length(unique(`Locus Tag`)), by=list(Source, GenomeID, SpeciesName)]
      setnames(gene_annot1, c("V1","Source"), c("CopyNum","KO"))
      all_gene_annots = rbind(all_gene_annots, gene_annot1)
    }
  }
  return(all_gene_annots)
}

####################### 3) Plotting Functions ##########################################
#Plotting functions

#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

#Modify axis breaks
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    signif_dig = 2
    #rounding_dig = ifelse(median(abs(x)) < 0.001|diff(range(x))< 0.005, 4, 2)
    #rounding_dig = ifelse(median(abs(x)) < 0.0005, 4, rounding_dig)
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    breaks = seq(signif(min(x)+d, signif_dig), signif(max(x)-d,signif_dig), length=n)
    breaks2 = signif(breaks, signif_dig)
    if(diff(range(breaks2)) < 0.006 & min(breaks2) > 0.01) breaks2 = signif(breaks, signif_dig+1)
    breaks2[abs(breaks2)< 0.0001] = 0
    if(all(breaks2 < 0)|all(breaks2 > 0)) return(breaks2) else return(sort(c(0, breaks2)))
  }
}

#Met fluxes plot
make_met_plot = function(metLab, sub_fluxes_plot, version = "subFluxes"){
  if(version=="subFluxes"){
    put_fluxes = ggplot(sub_fluxes_plot[niceLab==metLab & SpeciesName != "Inflow"], aes(x=SimRun, y = cumulFlux, col = SpeciesName, shape = type2))+ 
      geom_point(size=0.8)  + facet_wrap(~type2, scales = "free", nrow = 2) + theme_cowplot() +theme(strip.text.y = element_text(),axis.title.y = element_blank(), 
                                                                                                     axis.text.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(size=8), strip.background = element_blank(), 
                                                                                                     strip.text.x = element_blank(),  axis.ticks = element_blank(), axis.line = element_blank(), plot.background = element_blank(), 
                                                                                                     panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype="solid"), 
                                                                                                     legend.title = element_blank(), plot.margin = margin(0.15, 0.15, 0.15, 0.05, "inches")) + 
      scale_color_manual(values = col_spec3) + scale_x_continuous(expand = c(0.005,0)) + scale_y_continuous(breaks = equal_breaks()) + 
      guides(col = F, shape = F) 
  } else {
    put_fluxes = ggplot(sub_fluxes_plot[niceLab==metLab & SpeciesName != "Inflow"], aes(x=SimRun, y = cumulFlux, col = Code, shape = type2))+ geom_point(size=0.8)  + facet_wrap(~type2, scales = "free", nrow = 2) + theme_cowplot() +theme(strip.text.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(size=8), strip.background = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.background = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype="solid"), legend.title = element_blank(), legend.text = element_text(size=8), plot.margin = margin(0.15, 0.15, 0.15, 0.15, "inches")) + scale_color_manual(values = col_spec1) + scale_x_continuous(expand = c(0.005,0)) + scale_y_continuous() + scale_shape_manual(values = c(17,16)) + guides(col = F, shape = F)
  }
  return(put_fluxes) 
}

#Metabolite contributions plot
make_contribs_plot = function(metLab, met_sub_dat){
  put_contribs = ggplot(met_sub_dat[niceLab==metLab],  aes(y=VarShare, x = SpeciesName, fill = SpeciesName)) + geom_bar(stat = "identity") + 
    facet_wrap(~niceLab, scales = "free", nrow = 5) + scale_fill_manual(values = col_spec2) + geom_abline(intercept = 0, slope = 0, linetype = 2) + 
    theme(strip.background = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(), 
          legend.title = element_blank(), strip.text = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.15, 0.5, 0.05, 0.27, "inches")) + 
    ylab("Contribution to variance") + scale_y_continuous() + coord_flip() + guides(fill = F) + ylab("")#
  return(put_contribs) 
}


###################### 4) Miscellaneous #################################################
##########Calculate a correlation matrix (data.table)

basic_correlation_matrix = function(fake_spec, fake_mets, method="pearson"){
  all_dats = merge(fake_spec, fake_mets, by = "Sample", allow.cartesian = T) #big table
  cor_results = all_dats[,cor.test(value.x, value.y, method = method, use = "complete.obs")[c("estimate", "p.value")], by = list(Species, compound)]
  return(cor_results)
}

########### Take a classification table and calculate an ROC curve 
make_ROC = function(outcome_dat, outcome_var, measure_var){
  ROC_datAll = data.table()
  df = outcome_dat[,list(get(outcome_var), get(measure_var))]
  setnames(df, c("survived", "pred"))
  ROC_dat = data.table(calculate_roc(df[!is.na(survived)], 1,1, 5000))
  ROC_dat = unique(ROC_dat[,list(tpr, fpr)])
  return(ROC_dat) 
}

############ Get AUC from ROC curve
get_AUC = function(ROC_dat){
  return(sum(sapply(1:nrow(ROC_dat), function(x){ 
    return(ROC_dat[x,tpr]*(ROC_dat[x,fpr]-ROC_dat[x+1, fpr]))
  }), na.rm=T))
}

############## ROC calculation
calculate_roc <- function(df, cost_of_fp, cost_of_fn, n=100) {
  tpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == 1) / sum(df$survived == 1)
  }
  
  fpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == 0) / sum(df$survived == 0)
  }
  
  cost <- function(df, threshold, cost_of_fp, cost_of_fn) {
    sum(df$pred >= threshold & df$survived == 0) * cost_of_fp + 
      sum(df$pred < threshold & df$survived == 1) * cost_of_fn
  }
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn))
  
  return(roc)
}

##### Process FBA simulation output, fluxes is whether output_details was true for the simulations
processFBAOutput = function(matFile, fluxes = F, all_rxns_by_spec = NULL, agora = F){ 
  matObject = (readMat(matFile))
  if(fluxes) extraInfo = matObject$cocultures.detail
  
  matObject = matObject$cocultures
  nruns = dim(matObject)[3]
  all_results = list()
  for(k in 1:nruns){
    if(agora == F){
      spec_all = unlist(matObject[[(nruns-1)*22+1]]) #species names
      met_list = unlist(matObject[[(nruns-1)*22+2]]) #metabolite names
      #Species
      specAbunds = data.table(matObject[[(k-1)*22+15]])
      specAbunds[,Species:=spec_all]
      specAbunds = melt(specAbunds, id.var = "Species", variable.name = "TimePoint")  
      specAbunds[,TimePoint:=as.numeric(gsub("V","",TimePoint))]
      #Remove species not actually present
      good_spec = specAbunds[,length(value[value!=0]), by=Species][V1 > 0][,Species]
      specAbunds = specAbunds[Species %in% good_spec]
      #Metabolites
      metAbunds = data.table(matObject[[(k-1)*22+16]])
      metAbunds[,Metabolite:=met_list]
      metAbunds = melt(metAbunds, id.var = "Metabolite", variable.name = "TimePoint")
      metAbunds[,TimePoint:=as.numeric(gsub("V","", TimePoint))]
      #This is wrong
      #metStarting = data.table(Metabolite = met_list, TimePoint = rep(0), value = matObject[[(nruns-1)*22+17]][,1])
      #metAbunds = rbind(metAbunds, metStarting)
      good_mets = metAbunds[,length(value[value!=0]), by=Metabolite][V1 > 0][,Metabolite] #Nonzero at any time point
      metAbunds = metAbunds[Metabolite %in% good_mets]
      #Growth Rates
      growthRates = data.table(rbind(sapply(matObject[[(k-1)*22+18]], function(x){ return((x[[1]]))})))
      setnames(growthRates, spec_all)
      growthRates[,TimePoint:=1:nrow(growthRates)]
      growthRates = melt(growthRates, id.var = "TimePoint", variable.name = "Species")
      growthRates = growthRates[Species %in% good_spec]
    } else {
      info = matObject[,,k]
      spec_all = unlist(info$species)
      met_list = unlist(info$mets)
      #Species
      specAbunds = data.table(info$bio)
      specAbunds[,Species:=spec_all]
      specAbunds = melt(specAbunds, id.var = "Species", variable.name = "TimePoint")  
      specAbunds[,TimePoint:=as.numeric(gsub("V","",TimePoint))]
      #Remove species not actually present
      good_spec = specAbunds[,length(value[value!=0]), by=Species][V1 > 0][,Species]
      specAbunds = specAbunds[Species %in% good_spec]
      #Metabolites
      metAbunds = data.table(info$ctns)
      metAbunds[,Metabolite:=met_list]
      metAbunds = melt(metAbunds, id.var = "Metabolite", variable.name = "TimePoint")
      metAbunds[,TimePoint:=as.numeric(gsub("V","", TimePoint))]
      #Growth Rates
      growthRates = data.table(rbind(sapply(info$mu.model, function(x){ return((x[[1]]))})))
      setnames(growthRates, spec_all)
      growthRates[,TimePoint:=1:nrow(growthRates)]
      growthRates = melt(growthRates, id.var = "TimePoint", variable.name = "Species")
      growthRates = growthRates[Species %in% good_spec]
    }
    #Fluxes, need info on reactions for each species - now we have to fix this
    if(fluxes){
      if(is.null(all_rxns_by_spec)) stop("Model information for each species is required")
      #specRxnAbunds = FBA_picrust(specAbunds, )
      rxn_fluxes = rbindlist(lapply(1:length(spec_all), function(x){
        return(data.table(extraInfo[,,k]$allflux[[x]][[1]], Rxn = all_rxns_by_spec[[x]], Species = spec_all[x])) }) )
      rxn_fluxes = melt(rxn_fluxes, id.vars = c("Rxn", "Species"), variable.name = "TimePoint")
      rxn_fluxes[,TimePoint:=as.numeric(gsub("V","",TimePoint))+1]
      #Need species/rxn abundance to get cumulative flux
      #rxn_fluxes[,cumulFlux:=cumsum(value), by=list(Species, Rxn)]
      
      glpk_fluxes = rbindlist(lapply(1:length(spec_all), function(x){
        return(data.table(extraInfo[,,k]$glpkflux[[x]][[1]], Rxn = all_rxns_by_spec[[x]], Species = spec_all[x])) }) )
      glpk_fluxes = melt(glpk_fluxes, id.vars = c("Rxn", "Species"), variable.name = "TimePoint")
      glpk_fluxes[,TimePoint:=as.numeric(gsub("V","",TimePoint))+1]
      
      models2medium = rbindlist(lapply(1:length(spec_all), function(x){
        return(data.table(extraInfo[,,k]$v.model2medium[[x]][[1]], Metabolite = met_list, Species = spec_all[x])) }) )
      models2medium = melt(models2medium, id.vars = c("Metabolite", "Species"), variable.name = "TimePoint")
      models2medium[,TimePoint:=as.numeric(gsub("V","",TimePoint))+1]
      all_results[[k]] = list(species = specAbunds, metabolites = metAbunds, growthRates = growthRates, fluxes = rxn_fluxes, glpk_fluxes = glpk_fluxes, model2medium = models2medium)
    } else {
      all_results[[k]] = list(species = specAbunds, metabolites = metAbunds, growthRates = growthRates)
    }
  }
  return(all_results)
}


############# Get info from model file
getModelInfo = function(matFile){
  return(readMat(matFile)[[1]][,,1])
}

########## Get metabolite concentration changes
getMetaboliteDeltas = function(all_mets){
  all_mets2 = copy(all_mets)
  all_mets2[,TimePointPrev:=TimePoint-1]
  all_mets2[,Sample:=NULL]
  all_mets_diff = merge(all_mets, all_mets2, by.x = c("SimRun", "Metabolite", "medium", "TimePoint"), by.y =c("SimRun","Metabolite", "medium", "TimePointPrev"), all = T)
  all_mets_diff[,MetChange:=value.y-value.x]
  setnames(all_mets_diff, c("value.y", "value.x"), c("value", "valuePrev"))
  #Get rid of extra TimePoint column
  all_mets_diff[,TimePoint:=NULL]
  all_mets_diff = all_mets_diff[!is.na(TimePoint)]
  all_mets_diff = all_mets_diff[TimePoint > 0]
  return(all_mets_diff)
}

########### Get correct total cumulative fluxes
#Including inflow and outflow
#Or only including inflow and dispersing outflow - in this case it will only return for the last time point currently
getCorrectTotalFluxes = function(met_fluxes, mets, media, species, growthRates, d = 0.0472, t = 0.25, disperse_outflow = T, max_tp = 0){
  #Species & growth rates
  spec_info = merge(species, growthRates, by=c("SimRun", "Species", "TimePoint"), all = T)
  spec_info[,e_mu_t_minus1:=(exp(value.y*t)-1)]
  
  #Combine relevant info
  met_fluxes[,TimePoint2:=TimePoint-1]
  
  if(!("medium" %in% names(met_fluxes) & "Metabolite" %in% names(met_fluxes))){
    stop("Merge mod2media with dictionary to get both metabolite IDs")
  }
  mod2media_info = merge(met_fluxes, spec_info, by.x = c("SimRun", "Species", "TimePoint2"), by.y = c("SimRun", "Species", "TimePoint"), all = T)
  setnames(mod2media_info, c("value", "value.x", "value.y"), c("mod2media_v", "speciesAbund", "growthRate"))
  mod2media_info[,B_e_mu_t_minus1:=speciesAbund*e_mu_t_minus1]
  mod2media_info[,v_over_mu:=mod2media_v/growthRate]
  mod2media_info[,totFlux:=v_over_mu*B_e_mu_t_minus1]  
  
  if(disperse_outflow == F){
    #Set up metabolite concentrations
    mets = getMetaboliteDeltas(mets)
    mets[,Outflow:=-1*valuePrev*d*t] #Opposite sign of inflow
    #Merge with media info
    if(!("mediaAmt" %in% names(media))){
      setnames(media, c("V1", "V3"), c("medium", "mediaAmt"))
    }
    mets = merge(mets, media[,list(medium, mediaAmt)], by="medium", all.x=T, all.y=F)
    mets[TimePoint==1, Inflow:=mediaAmt]
    mets[TimePoint > 1,Inflow:=mediaAmt*t]
    mets[is.na(Inflow), Inflow:=0]
    
    #Combine everything
    external_fluxes = melt(mets[,list(SimRun, TimePoint, Metabolite, medium, Inflow, Outflow)], id.vars = c("SimRun", "TimePoint", "Metabolite", "medium"), variable.name = "Species", value.name = "totFlux")
    internal_fluxes = mod2media_info[,list(SimRun, Species, TimePoint, Metabolite, medium, totFlux)]
    all_met_fluxes = rbind(external_fluxes, internal_fluxes, fill = T)
    
    #Get cumulative fluxes
    all_met_fluxes = all_met_fluxes[!is.na(Metabolite) & !is.na(totFlux) & !is.na(TimePoint)]
    all_met_fluxes[order(TimePoint), cumulFlux:=cumsum(totFlux), by=list(SimRun, Species, Metabolite)]
    return(all_met_fluxes)
  } else {
    if(max_tp == 0) max_tp = met_fluxes[,max(TimePoint)]
    tot1 = mod2media_info[TimePoint <= max_tp,list(SimRun, Species, TimePoint, Metabolite, medium, totFlux)]
    tot1[,adjFlux:=totFlux*(1-d*t)^(max_tp - TimePoint)]
    media_adj = data.table(Species = "Inflow", unique(mod2media_info[TimePoint <= max_tp,list(medium, Metabolite, SimRun,TimePoint)]))
    media_adj = rbind(media_adj, data.table(TimePoint=1, Species = "Inflow", unique(mod2media[,list(SimRun,Metabolite, medium)])), fill = T)
    media_adj = merge(media_adj, media[,list(V1, V3)], by.x = "medium", by.y = "V1", all.x = T)
    media_adj[is.na(V3), V3:=0]
    media_adj[,totFlux:=ifelse(TimePoint==1, V3, V3*t)]
    media_adj[,adjFlux:=totFlux*(1-d*t)^(max_tp - TimePoint)]
    all_met_fluxes = rbind(tot1, media_adj, fill = T)
    all_met_fluxes[order(TimePoint), cumulFlux:=cumsum(adjFlux), by=list(SimRun, Species, Metabolite)]
    all_fluxes_final= all_met_fluxes[TimePoint==max_tp]
    return(all_fluxes_final)
  }
}

########### Get correct total cumulative fluxes
#No accounting for outflow - just getting total activity
getCorrectTotalRxnFluxes = function(fluxes, species, growthRates, d = 0.0472, t = 0.25){
  #Species & growth rates
  spec_info = merge(species, growthRates, by=c("SimRun", "Species", "TimePoint"), all = T)
  spec_info[,e_mu_t_minus1:=(exp(value.y*t)-1)]
  
  #Combine relevant info
  fluxes[,TimePoint2:=TimePoint-1]
  
  fluxes_info = merge(fluxes, spec_info, by.x = c("SimRun", "Species", "TimePoint2"), by.y = c("SimRun", "Species", "TimePoint"))#, all = T)
  setnames(fluxes_info, c("value", "value.x", "value.y"), c("mod2media_v", "speciesAbund", "growthRate"))
  fluxes_info[,B_e_mu_t_minus1:=speciesAbund*e_mu_t_minus1]
  fluxes_info[,v_over_mu:=mod2media_v/growthRate]
  fluxes_info[,totFlux:=v_over_mu*B_e_mu_t_minus1]  
  
  all_rxn_fluxes = fluxes_info[!is.na(Rxn) & !is.na(totFlux) & !is.na(TimePoint)]
  all_rxn_fluxes[order(TimePoint), cumulFlux:=cumsum(totFlux), by=list(SimRun, Species, Rxn)]
  return(all_rxn_fluxes)
}

########## Get species table
make_spec_codes = function(){
  spec_codes = data.table(Species = c("ErectaleAGORA", "CaerofaciensAGORA", "DpigerAGORA", "BhydrogenotrophicaAGORA", "CsymbiosumAGORA", "MformatexigensAGORA", "EcoliAGORA", "BovatusAGORA", "BthetaiotaomicronAGORA", "BcaccaeAGORA"), Code = c("Er", "Ca", "Dp", "Bh", "Cs", "Mf", "Ec", "Bo", "Bt", "Bc"), SpeciesName = c("Eubacterium rectale", "Collinsella aerofaciens", "Desulfovibrio piger", "Blautia hydrogenotrophica", "Clostridium symbiosum", "Marvinbryantia formatexigens", "Escherichia coli", "Bacteroides ovatus", "Bacteroides thetaiotaomicron", "Bacteroides caccae"), Class = c("Clostridia", "Actinobacteria", "Deltaproteobacteria", "Clostridia", "Clostridia","Clostridia", "Gammaproteobacteria", "Bacteroidia", "Bacteroidia", "Bacteroidia"), Phylum = c("Firmicutes", "Actinobacteria", "Proteobacteria", "Firmicutes", "Firmicutes", "Firmicutes", "Proteobacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidetes"), Description = c("Eubacterium rectale ATCC 33656", "Collinsella aerofaciens ATCC 25986", "Desulfovibrio piger ATCC 29098", "Blautia hydrogenotrophica DSM 10507", "Clostridium symbiosum WAL-14163", "Marvinbryantia formatexigens I-52 DSM 14469", "Escherichia coli K-12 MG1655", "Bacteroides ovatus ATCC 8483", "Bacteroides thetaiotaomicron VPI-5482", "Bacteroides caccae ATCC 43185"))
  spec_codes = rbind(spec_codes, data.table(Species = c("Inflow", "Outflow"), Code = c("in", "out")), fill = T)
  spec_codes[is.na(SpeciesName), SpeciesName:=Species]
  return(spec_codes)
}

############### Fix met names to be clean/compatible - AGORA dictionary
fix_met_names = function(met_data, col_name = "Metabolite"){
  all_mets[Metabolite=='L-lysinium(1+)', Primary:="L-lysine"]
  dictionary[,Primary:=gsub(" \\(.*\\)$","", Primary)]
  dictionary[,Primary:=gsub("\\(.*\\)$","", Primary)]
  dictionary[Primary=='L-argininium', Primary:="L-arginine"]
  dictionary[Primary=="proton", Primary:="H+"]
  dictionary[Primary=="hydrogenphosphate", Primary:="Orthophosphate"]
  dictionary[grepl("Amylopectin", Primary), Primary:="Amylopectin"]
  setnames(dictionary, "Primary", "Metabolite")
  
}

####################### Build dataset with processed reaction IDs from multiple sets of simulations, must all share the same file name structure 
process_FBA_runs = function(mod_path, results_path, file_id, dictionary, spec_codes, out_prefix){
  
  #Load models
  all_mods = list()
  for(j in 1:length(spec_codes[,Species])){
    all_mods[[j]] = getModelInfo(paste0(mod_path, spec_codes[j,Species], ".mat"))
  }
  
  all_rxns_by_spec = lapply(all_mods, function(x){ return(unlist(x$rxns))}) #Includes biomass0, does not include fixed conflict rxns
  
  mimosa_pieces = build_model_components(all_mods, remove_rev = T, missing_rxns = T, spec_codes, agora = T)
  
  all_spec = paste(spec_codes[,Code], collapse="")
  ten_spec_vary = list()
  for(j in 1:length(input_files)){
    #Process in sequential order
    spec_vary = processFBAOutput(paste0(results_path, "/", file_id, "_", j, ".mat"), fluxes = T, all_rxns_by_spec, agora = T)
    ten_spec_vary = append(ten_spec_vary, spec_vary)
  }
  
  ten_spec_vary = lapply(ten_spec_vary, fix_rxn_ids, specRxns_melt_rev = mimosa_pieces[[9]])
  ten_spec_vary = lapply(ten_spec_vary, get_total_and_cumul_fluxes, specRxns = mimosa_pieces[[8]], spec_codes)
  #Split rev rxn IDs
  fixed_rxns = fix_rev_rxns(ten_spec_vary, mimosa_pieces)
  ten_spec_vary_revFix = fixed_rxns[[1]]
  mimosa_pieces_revFix = fixed_rxns[[2]]
  all_mets = rbindlist(lapply(1:length(ten_spec_vary_revFix), function(x){
    met_sub = ten_spec_vary_revFix[[x]]$metabolites
    met_sub[,SimRun:=x]
    return(met_sub)
  }))
  all_mets = merge(all_mets, dictionary[,list(Primary, medium)], by.x="Metabolite", by.y = "Primary", all.x=T, all.y = F)
  
  all_species = rbindlist(lapply(1:length(ten_spec_vary_revFix), function(x){foo = ten_spec_vary_revFix[[x]]$species 
  foo[,SimRun:=x]
  return(foo)}))
  
  growthRates = rbindlist(lapply(1:length(ten_spec_vary_revFix), function(x){
    foo = ten_spec_vary_revFix[[x]]$growthRates
    foo[,SimRun:=x]
    return(foo)
  }))
  
  mod2media = rbindlist(lapply(1:length(ten_spec_vary_revFix), function(x){
    foo = ten_spec_vary_revFix[[x]]$model2medium
    foo[,SimRun:=x]
    return(foo)
  }))
  
  save(all_species, file = paste0(out_prefix, "_allSpecies.rda"))
  save(mimosa_pieces_revFix, file = paste0(out_prefix, "_mimosa_pieces.rda"))
  #save(ten_spec_vary_revFix, file = paste0(out_prefix, "_resultsList.rda"))
  save(all_fluxes, file = paste0(out_prefix, "_allFluxes.rda"))
  save(all_mets, file = paste0(out_prefix, "_allMets.rda"))
  save(growthRates, file = paste0(out_prefix, "_growthRates.rda"))
  save(mod2media, file = paste0(out_prefix, "_mod2media.rda"))
}

########### Add species reaction abundance to flux table
get_specRxn_flux_data = function(run_data, specRxns, spec_codes){
  spec_data = run_data$species
  spec_data = dcast.data.table(spec_data, Species~TimePoint, value.var = "value", fun.aggregate = sum)
  #Get gene abundances
  gene_data_by_spec = rbindlist(lapply(1:nrow(spec_data), function(x){
    spec_gene_data = melt(FBA_picrust(spec_data[Species==spec_codes[x,Species]], specRxns, spec_codes), id.var = "KO", variable.name = "TimePoint", value.name = "specRxnAbund")
    spec_gene_data[,Species:=spec_codes[x,Species]]
    return(spec_gene_data[specRxnAbund !=0])
  }))
  setnames(gene_data_by_spec, "KO", "Rxn")
  gene_data_by_spec[,TimePoint:=as.numeric(TimePoint)]
  run_data$fluxes = merge(run_data$fluxes, gene_data_by_spec, by = c("Species","Rxn","TimePoint"), all = T)
  run_data$fluxes = run_data$fluxes[Rxn !="biomass0"& !grepl("EX_", Rxn)] # 
  run_data$fluxes[TimePoint==1,value:=0]
  run_data$fluxes[,totalFlux:=value*specRxnAbund]
  run_data$fluxes[order(TimePoint),cumulFlux:=cumsum(totalFlux),by=list(Species,Rxn)] #This is a placeholder (calculated correctly later)
  return(run_data)
}


############### Convert community stoichiometric matrix (e.g. MIMOSA) to edge list
emm_to_edge_list = function(emm){ 
  net_melted = melt(emm, id.var = "Compound")[value != 0]
  net_melted[,Prod:=ifelse(value > 0, Compound,0)]
  net_melted[,Reac:=ifelse(value < 0, Compound,0)]
  all_rxn_ids = net_melted[,unique(as.character(variable))]
  edge_list = data.table()
  
  for(k in 1:length(all_rxn_ids)){
    rxn_sub = net_melted[variable==all_rxn_ids[k]]
    if(nrow(rxn_sub[Prod !=0 ]) > 0 & nrow(rxn_sub[Reac != 0]) > 0){
      edge_list_sub = data.table(expand.grid(rxn_sub[,unique(Reac[Reac !=0])], rxn_sub[,unique(Prod[Prod != 0])]))
      setnames(edge_list_sub, c("Reac", "Prod"))
    } else {
      edge_list_sub = rxn_sub[,list(Reac, Prod)]
      edge_list_sub[Reac==0, Reac:=NA]
      edge_list_sub[Prod==0, Prod:=NA]
    }
    edge_list_sub[,KO:=all_rxn_ids[k]]
    edge_list_sub[,stoichReac:=sapply(Reac, function(x){ return(rxn_sub[Compound==x,abs(value)])})]
    edge_list_sub[,stoichProd:=sapply(Prod, function(x){ return(rxn_sub[Compound==x,value])})]
    edge_list = rbind(edge_list, edge_list_sub)
  }
  return(edge_list)
}

#Get network distance between 2 compounds (does not account for directionality)
get_net_dist = function(c1, c2, allnet, currency_mets = data.table()){
  if(c1==c2) return(0)
  match_rxns = allnet[Prod==c1 | Reac==c1]
  netdist = 1
  while(!(c2 %in% match_rxns[,Prod]) & !(c2 %in% match_rxns[,Reac]) & netdist < 20){
    netdist = netdist + 1
    if(!(c1 %in% currency_mets) & !(c2 %in% currency_mets)){
      all_comps = unique(c(match_rxns[ProdC==F,Prod], match_rxns[ReacC==F,Reac]))
      match_rxns = allnet[Prod %in% all_comps | Reac %in% all_comps]
    } else {
      all_comps = unique(c(match_rxns[,Prod], match_rxns[,Reac]))
      match_rxns = allnet[Prod %in% all_comps | Reac %in% all_comps]
    }
  }
  return(netdist)
}

################ Make media files with noisy quantities, write to same directory as original file
make_media_vary = function(orig_media_file, data_dir, used_mets, nsamps = 61, sd_constants = c(0.5, 1, 2, 3, 4, 5, 8, 10)){
  main_media = fread(paste0(data_dir, orig_media_file))
  media_mets = main_media[,V1]
  if(any(grepl("[e]", used_mets, fixed = T))){
    used_mets = used_mets[gsub("[env]", "[e]", Metabolite, fixed = T)]
  }
  media_mets = media_mets[media_mets %in% used_mets]
  
  for(k in 1:length(sd_constants)){
    for(j in 1:nsamps){
      media2 = copy(main_media)
      media2[, V3:=sapply(V3, function(x){ return(rnorm(1, mean = x, sd = x*sd_constants[k]/100))})]
      write.table(media2, file = paste0(data_dir, "FaithMediaAGORA_", sd_constants[k], "_", j, ".csv"), quote=F, row.names = F, col.names = F, sep = ",")
    }
  }
}

################### Generate metabolite categories
make_path_key = function(met_key){
  load("../MIMOSA/mimosa/R/sysdata.rda")
  met_key = fread(paste0(homedir, "DATA/hmdb_metabolite_classes.txt"))
  met_key = met_key[Name != "None"]
  met_key[,altName:=gsub("ic acid$", "ate",Name)]
  met_key[,Name:=tolower(Name)]
  met_key[,altName:=tolower(altName)]
  
  bad_mets = dictionary_sub[!(tolower(Metabolite) %in% met_key[,tolower(Name)]) & !(tolower(Metabolite)) %in% met_key[,tolower(altName)], Metabolite]
  
  bad_mets_kegg = lapply(bad_mets, keggFind, database="compound")
  bad_mets_kegg = lapply(bad_mets_kegg, function(x){
    return(sapply(x, function(y){ c(strsplit(y, split = "; "))}))
  })
  bad_mets_match = lapply(1:length(bad_mets_kegg), function(x){
    if(length(bad_mets_kegg[[x]]) !=0){
      test1 = sapply(1:length(bad_mets_kegg[[x]]), function(y){
        if(any(tolower(bad_mets_kegg[[x]][[y]])==tolower(bad_mets[x]))) return(names(bad_mets_kegg[[x]][y])) else return(NA)
      })
      if(any(!is.na(test1))) return(test1[!is.na(test1)]) else return(NA)
    } else return(NA)
  })
  bad_mets_match = gsub("cpd:", "", unlist(bad_mets_match))
  bad_mets = data.table(Metabolite = bad_mets, KEGG = bad_mets_match)
  bad_mets[Metabolite=="octadecanoate", KEGG:="C01530"]
  bad_mets[Metabolite=="Zinc", KEGG:="C00038"]
  bad_mets[grepl("Isobutyrate", Metabolite), KEGG:="C02632"]
  bad_mets[grepl("Menaquinone", Metabolite), KEGG:="C00828"]
  bad_mets[grepl("Diaminoheptanedioate", Metabolite), KEGG:="C00680"]
  bad_mets = merge(bad_mets, met_key, by = "KEGG", all.x=T, all.y=F)
  
  path_key = met_key[Name %in% dictionary_sub[,tolower(Metabolite)]]
  path_key_alt = met_key[altName %in% dictionary_sub[,tolower(Metabolite)]]
  setnames(path_key_alt, c("Name", "altName"), c("altName", "Name"))
  path_key = unique(rbind(path_key, path_key_alt, fill = T))
  path_key = unique(rbind(path_key, bad_mets, fill = T))
  path_key[is.na(SubClass), SubClass:="Other"]
  path_key[is.na(Class), Class:="Other"]
  path_key[Class=="None", Class:="Other"]
  path_key[SubClass=="None", SubClass:='Other']
  path_key[is.na(Name), Name:=tolower(Metabolite)]
  path_key[,niceLab:=ifelse(!is.na(Metabolite), tolower(Metabolite), tolower(Name))]
  path_key_bad = path_key[,which(duplicated(niceLab))]
  path_key = path_key[-path_key_bad]
  return(path_key)
}



