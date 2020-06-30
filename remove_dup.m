[Data,Dir] = uigetfile('.txt');
FPdata = load(strcat(Dir,Data));
intensity_dup = FPdata(:,1).';
[lambda, idx] = unique(FPdata(:,2).');
lambda = lambda.';
error_dup = FPdata(:,3).';

intensity = zeros(size(idx));
error = zeros(size(idx));
for n = 1:length(idx)
    intensity(n) = intensity_dup(idx(n));
    error(n) = error_dup(idx(n));
end

M = [intensity lambda error];

dlmwrite('FP 5.32mm WG high res 1pm (1).txt',M,'delimiter','\t','precision',10)


