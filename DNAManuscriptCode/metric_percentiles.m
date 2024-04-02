seq_list_string=strings;

%Use this to get median, top 0.1%, and bottom 0.1%
%Also prints median sequence
for i=1:10000
    operators = {'A', 'T', 'G', 'C'};
    operators = {'A', 'T', 'G', 'C', 'A','T'};
    operators = {'A', 'T', 'G', 'C', 'G','C'};
    indexes = randi(length(operators), 1, 100);
    seq=cell2mat(operators(indexes));
    seq_list_string(i)=seq;
    %Choose index of metric of interest
    %Call runningforjavacgDNA for cgDNA rather than crystal statistics
    list_results(i,1)=runningforjava(seq,1); 
end
list_results=sortrows(list_results);
h=size(list_results);
actualsize=h(1,1);
list_results(actualsize*(500/1000),1)
list_results(actualsize*(999/1000),1)
list_results(actualsize*(1/1000),1)
seq_list_string(list_results(actualsize*(500/1000),2))
