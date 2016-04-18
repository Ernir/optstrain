function verbosePrint(string, verbosity)
% If verbosity is true, very prominently display the string.
if verbosity
    strLength = length(string);
    disp(char(ones(1,strLength + 2)*double('%')))
    disp(['% ' string])
    disp(char(ones(1,strLength + 2)*double('%')))
end
end