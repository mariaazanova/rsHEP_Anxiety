function [n_ext, hr] = get_cers_extrasystoles(cers_file, subjid)
% GET_CERS_EXTRASYSTOLES  V + S beat count and heart rate (beats/min) of one participant in the CER-S
%   summary table (row ECG_<subjid>_Ruhe); both NaN without a row or CER-S result. '-' means none.
%   hr = # QRS / analyzed time.
lines = splitlines(fileread(cers_file));
hdr = strtrim(strsplit(lines{1}, ';'));
names = cellfun(@(l) strtrim(strtok(l, ';')), lines, 'UniformOutput', false);
row = find(strcmpi(names, ['ECG_', subjid, '_Ruhe.hea.csv']));
n_ext = NaN; hr = NaN;
if isempty(row), return, end
f = strtrim(strsplit(lines{row}, ';'));
val = @(name) f{strcmp(hdr, name)};
if numel(f) < numel(hdr) || isnan(str2double(val('# QRS'))), return, end
n_ext = sum(str2double(strrep({val('V Beats'), val('S Beats')}, '-', '0')));
t = str2double(regexp(val('Analyzed Time'), '^(\d+):(\d+):(\d+)', 'tokens', 'once'));
hr = str2double(val('# QRS')) / ([60 1 1/60] * t(:));
end
