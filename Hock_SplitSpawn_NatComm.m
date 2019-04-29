% Code to calculate the connectivity metrics for the article
% "Split spawning increases robustness of coral larval supply and
% inter-reef connectivity", by Hock et al. (2019) Nature Communications
% Author & copyright: Karlo Hock, University of Queensland. 2018

% The matrices model yearly/seasonal inter-reef connectivity of Acropora corals on the Great Barrier
% Reef over 7 seasons and for 9 variants of latitidunal sectors; see the article for details
% The code needs to load three data files: ACROPsplits.mat, ACROPm1.mat,
% ACROPm2.mat; and access two function files: scalebysize, zerodiag

% Data files can be downloaded from 10.5281/zenodo.2653244

% Split spawning%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, calculate everything for split spawning scenarios, using 60:40
% ratio for split spawning on the Great Barrier Reef from Willis et al 1985

acrowithsplits=struct('srcs',[],'sply',[]);
load ACROPsplits.mat;
for sector=1:9% calculate everything for different sector boundaries to use later in sensitivity analysis
    n_sources=zeros(3806,7);% number of sources per year
    tot_supply=zeros(3806,7);% total external supply per year
    for yr=1:7
        thisnet=zerodiag(cell2mat(ACROS(sector,yr)));% remove diagonal so that only external supply is considered
        thisnet=scalebysize(thisnet,sizes);% scale supply according to the size of the source reef polygon
        for j=1:3806
            n_sources(j,yr)=nnz(thisnet(:,j));
            tot_supply(j,yr)=sum(thisnet(:,j));
        end
    end
    yearly_totsply=zeros(3806,7);
    for reef=1:3806% normalise supply to the largest value of that year 
        for year=1:7
            yearly_totsply(reef,year)=tot_supply(reef,year)/max(tot_supply(:,year));
        end
    end
    tot_supply(:,1:7)=yearly_totsply;
    for rf=1:3806% calculate per reef metrics
        tot_supply(rf,8)=mean(tot_supply(rf,1:7));% fig 2A
        tot_supply(rf,9)=std(tot_supply(rf,1:7))/mean(tot_supply(rf,1:7));% fig 2B
        n_sources(rf,8)=mean(n_sources(rf,1:7));% fig 3A
        n_sources(rf,9)=std(n_sources(rf,1:7))/mean(n_sources(rf,1:7));% fig 3B
        n_sources(rf,10)=7-nnz(n_sources(rf,1:7));% fig 4A
        consecyrs=n_sources(rf,1:7);
        n_sources(rf,11)=max(diff([0 (find(~(consecyrs<1))) numel(consecyrs)+1]) -1);% for fig 5A later
    end
    
    acrowithsplits(sector).srcs=n_sources;
    acrowithsplits(sector).sply=tot_supply;

end

fig2A=acrowithsplits(5).sply(:,8);
fig2B=acrowithsplits(5).sply(:,9);
fig3A=acrowithsplits(5).srcs(:,8);
fig3B=acrowithsplits(5).srcs(:,9);
fig4A=acrowithsplits(5).srcs(:,10);

clear ACROS;% clear memory so as to not choke an average machine


% First month%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Then, calculate everything as if coral spawned only in the first month; same
% procedures as with the split spawning case, but now using the matrices only for the
% first spawning month

acromonth1=struct('srcs',[],'sply',[]);
load ACROPm1.mat;
for sector=1:9
    n_sources=zeros(3806,7);
    tot_supply=zeros(3806,7);
    for yr=1:7
        thisnet=zerodiag(cell2mat(ACRO1(sector,yr)));
        thisnet=scalebysize(thisnet,sizes);
        for j=1:3806
            n_sources(j,yr)=nnz(thisnet(:,j));
            tot_supply(j,yr)=sum(thisnet(:,j));
        end
    end
    yearly_totsply=zeros(3806,7);
    for reef=1:3806
        for year=1:7
            yearly_totsply(reef,year)=tot_supply(reef,year)/max(tot_supply(:,year));
        end
    end
    tot_supply(:,1:7)=yearly_totsply;
    for rf=1:3806
        tot_supply(rf,8)=mean(tot_supply(rf,1:7));
        tot_supply(rf,9)=std(tot_supply(rf,1:7))/mean(tot_supply(rf,1:7));
        n_sources(rf,8)=mean(n_sources(rf,1:7));
        n_sources(rf,9)=std(n_sources(rf,1:7))/mean(n_sources(rf,1:7));
        n_sources(rf,10)=7-nnz(n_sources(rf,1:7));%fig 4B
        consecyrs=n_sources(rf,1:7);
        n_sources(rf,11)=max(diff([0 (find(~(consecyrs<1))) numel(consecyrs)+1]) -1);%for fig 5A later
    end
    
    acromonth1(sector).srcs=n_sources;
    acromonth1(sector).sply=tot_supply;

end

fig4B=acromonth1(5).srcs(:,10);% this shows years without external supply when only the first month of spawning is taken into account

clear ACRO1;% clear memory so as to not choke an average machine


% Second month%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finally, calculate everything as if coral spawned only in the second month; same
% procedures as with the split spawning case, but now using the matrices only for the
% second spawning month

acromonth2=struct('srcs',[],'sply',[]);
load ACROPm2.mat;
for sector=1:9
    n_sources=zeros(3806,7);
    tot_supply=zeros(3806,7);
    for yr=1:7
        thisnet=zerodiag(cell2mat(ACRO2(sector,yr)));
        thisnet=scalebysize(thisnet,sizes);
        for j=1:3806
            n_sources(j,yr)=nnz(thisnet(:,j));
            tot_supply(j,yr)=sum(thisnet(:,j));
        end
    end
    yearly_totsply=zeros(3806,7);
    for reef=1:3806
        for year=1:7
            yearly_totsply(reef,year)=tot_supply(reef,year)/max(tot_supply(:,year));
        end
    end
    tot_supply(:,1:7)=yearly_totsply;
    for rf=1:3806
        tot_supply(rf,8)=mean(tot_supply(rf,1:7));
        tot_supply(rf,9)=std(tot_supply(rf,1:7))/mean(tot_supply(rf,1:7));
        n_sources(rf,8)=mean(n_sources(rf,1:7));
        n_sources(rf,9)=std(n_sources(rf,1:7))/mean(n_sources(rf,1:7));
        n_sources(rf,10)=7-nnz(n_sources(rf,1:7));
        consecyrs=n_sources(rf,1:7);
        n_sources(rf,11)=max(diff([0 (find(~(consecyrs<1))) numel(consecyrs)+1]) -1);%for fig 5A later
    end
    
    acromonth2(sector).srcs=n_sources;
    acromonth2(sector).sply=tot_supply;

end
clear ACRO2;% clear memory so as to not choke an average machine

% Ratios between months%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate ratios of metrics between the first and the second month to determine whether
% the metric would be skewed if spawning occured only once per year

ratiosply=zeros(3806,7);% ratio of external supply between two months per season
for reef=1:3806
    for year=1:7
        ratiosply(reef,year)=(0.4*acromonth2(5).sply(reef,year))/((0.6*acromonth1(5).sply(reef,year))+(0.4*acromonth2(5).sply(reef,year)));
    end
end
ratiosply(isnan(ratiosply))=0;
for sector=1:3806
    ratiosply(sector,8)=mean(ratiosply(sector,1:7));%fig2C
    ratiosply(sector,9)=std(ratiosply(sector,1:7))/mean(ratiosply(sector,1:7));
end

ratiosrcs=zeros(3806,7);% ratio of sources between two months per season
for reef=1:3806
    for year=1:7
        ratiosrcs(reef,year)=acromonth2(5).srcs(reef,year)/acrowithsplits(5).srcs(reef,year);
    end
end
ratiosrcs(isnan(ratiosrcs))=0;
for rf=1:3806
    ratiosrcs(rf,8)=mean(ratiosrcs(rf,1:7));%fig3C
    ratiosrcs(rf,9)=std(ratiosrcs(rf,1:7))/mean(ratiosrcs(rf,1:7));
end

fig2C=ratiosply(:,8);
fig3C=ratiosrcs(:,8);

% Difference between split spawning and a single (first) month in terms of
% years without supply and consecutive years without supply
difffails=acromonth1(5).srcs(:,10)-acrowithsplits(5).srcs(:,10);%fig4C
diffcfails=acromonth1(5).srcs(:,11)-acrowithsplits(5).srcs(:,11);%fig5A

fig4C=difffails;% total years without external supply
fig5A=diffcfails;% consecutive years without external supply


% Sensitivity analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate how much the metrics will vary if latitiudinal sectors are defined
% differently with alternative placement of spawning sector boundaries; 
% this is performed using the split spawning scenario
% see Supplementary Figure 1 of the article for sector boundaries

meansrcs=zeros(3806,9);% number of sources
meansply=zeros(3806,9);% amount of supply
meanfail=zeros(3806,9);% total years without external supply
meancfail=zeros(3806,9);% consecutive years without external supply
for sector=1:length(acrowithsplits)
    meansrcs(:,sector)=acrowithsplits(sector).srcs(:,8);
    meansply(:,sector)=acrowithsplits(sector).sply(:,8);
    meanfail(:,sector)=acrowithsplits(sector).srcs(:,10);
    meancfail(:,sector)=acrowithsplits(sector).srcs(:,11);
end

for rf=1:3806% run per reef calculations for each sector variant
    if any(meansrcs(rf,1:9)==0)
        if all(meansrcs(rf,1:9)==0)
            meansrcs(rf,10)=1;
        else
            meansrcs(rf,10)=0;
        end
    else
        meansrcs(rf,10)=min(meansrcs(rf,1:9))/max(meansrcs(rf,1:9));
    end
    if any(meansply(rf,1:9)==0)
        if all(meansply(rf,1:9)==0)
            meansply(rf,10)=1;
        else
            meansply(rf,10)=0;
        end
    else
        meansply(rf,10)=min(meansply(rf,1:9))/max(meansply(rf,1:9));
    end
    if any(meanfail(rf,1:9)==0)
        if all(meanfail(rf,1:9)==0)
            meanfail(rf,10)=1;
        else
            meanfail(rf,10)=0;
        end
    else
        meanfail(rf,10)=min(meanfail(rf,1:9))/max(meanfail(rf,1:9));
    end
    if any(meancfail(rf,1:9)==0)
        if all(meancfail(rf,1:9)==0)
            meancfail(rf,10)=1;
        else
            meancfail(rf,10)=0;
        end
    else
        meancfail(rf,10)=min(meancfail(rf,1:9))/max(meancfail(rf,1:9));
    end
end

% Sensitivity is calculated as the ratio of minimum and maximum values for a metric
% this gives the percentage difference that would result in a metric with
% alternative boundaries for latitudinal sectors; 
% sens1-sens4 values per metric are provided in Supplementary Table 2 of the article
mm1=mean(meansply(:,1:9),1);
sens1=max(mm1)/min(mm1);
mm2=mean(meansrcs(:,1:9),1);
sens2=max(mm2)/min(mm2);
mm3=mean(meanfail(:,1:9),1);
sens3=max(mm3)/min(mm3);
mm4=mean(meancfail(:,1:9),1);
sens4=max(mm4)/min(mm4);


% Dividing the GBR into sectors using variable sector boundaries%%%%%%%%%%%

% This code needs to load reef centroids from the centroids.mat file

% Centroids have been derived from GBRMPA GIS file available at 
% http://www.gbrmpa.gov.au/resources-and-publications/spatial-data-information-services

% The main analysis does not need to run this code as the sectors have been 
% predefined for it; nevertheless, this code shows the details of how the 
% sectors have been defined in the article and for the sensitvity analysis

% Defining the sector boundaries; see Supplementary Figure 1 of the article 
% for the visualisation of different sector boundaries
Nlines=[146.51 -14.028 143.082 -17.456; 147.375 -14.814 143.918 -18.275; 148.249 -15.65 144.821 -19.078];
Slines=[150.941 -17.707 147.475 -21.168; 151.743 -18.593 148.299 -22.037; 152.563 -19.429 149.118 -22.873];
nedge=[146.5 -7 137 -10];
sedge=[152 -28 160 -20];
%1 - north, 2 - central, 3 - south

load centroids.mat;
sectors4sens=zeros(3806,9);% Container to hold different subdivisions of the GBR's 3806 reefs
cntr=1;
for i=1:size(Nlines,1)
    for j=1:size(Slines,1)
        [npolyX,npolyY]=poly2cw([nedge(1,1) Nlines(i,1) Nlines(i,3) nedge(1,3)],[nedge(1,2) Nlines(i,2) Nlines(i,4) nedge(1,4)]);
        npoly=[transpose(npolyX) transpose(npolyY)];
        [cpolyX,cpolyY]=poly2cw([Nlines(i,1) Slines(j,1) Slines(j,3) Nlines(i,3)],[Nlines(i,2) Slines(j,2) Slines(j,4) Nlines(i,4) ]);
        cpoly=[transpose(cpolyX) transpose(cpolyY)];
        [spolyX,spolyY]=poly2cw([Slines(j,3) sedge(1,1) sedge(1,3) Slines(j,1)],[Slines(j,4) sedge(1,2) sedge(1,4) Slines(j,2)]);
        spoly=[transpose(spolyX) transpose(spolyY)];
        for r=1:3806
            if inpoly(centr(r,:),npoly)
                sectors4sens(r,cntr)=1;
            end
            if inpoly(centr(r,:),cpoly)
                sectors4sens(r,cntr)=2;
            end
            if inpoly(centr(r,:),spoly)
                sectors4sens(r,cntr)=3;
            end
        end
        cntr=cntr+1;
    end
end
