function data = readJson(filename)
fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);

data = JSON.parse(str);
end