function outEEG =  parload(fname)
  x = load(fname)
  outEEG = x.x;
  
end