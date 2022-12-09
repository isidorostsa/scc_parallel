%% Data containers

% datasets: celegansneural, 
% datasets: celegansneural, foldoc, language, eu-2005, wiki-topcats, sx-stackoverflow, wb-edu, indochina-2004, uk-2002, arabic-2005, uk-2005, twitter7
% we want to access the execution time of each of those files, we can load it,
% create the map
% and then access it


implementations = ["serial", "OpenMP", "OpenCilk", "pthread"]

datasets = ["celegansneural", "foldoc", "language", "eu-2005", "wiki-topcats", "sx-stackoverflow", "wb-edu", "indochina-2004", "uk-2002", "arabic-2005", "uk-2005", "twitter7"];

nthreads = [1, 2, 4, 6, 8, 10, 12]

% each implementation has a map of files and their execution times for each number of threads
% we can access it as follows (for example):
% impl_times = execution_times("OpenMP");
% data_impl_times = impl_times("celegansneural")
% data_impl_td_times = data_implt_times(1)
% returns the execution time of celegansneural with 1 thread and OpenMP

% create a database of execution times
% columns: implementations (OpenMP, OpenCilk, pthread) 
%           files (celegansneural, foldoc, language, eu-2005, wiki-topcats, sx-stackoverflow, wb-edu, indochina-2004, uk-2002, arabic-2005, uk-2005, twitter7)
%           nthreads (1, 2, 4, 6, 8, 10, 12)

