function prefix = get_system_prefix(ELEM_TYPE)

[s,w] = system('hostname');
assert(s == 0);
  
if (strcmp(w, sprintf('ozzy.csl.uiuc.edu\n')) || ...
    strcmp(w, sprintf('samurai.local\n')))
  prefix = 'DYLD_LIBRARY_PATH=/usr/local/lib';
elseif(strcmp(w, sprintf('wutang\n')) || ...
       strcmp(w, sprintf('bacchus\n')) || ...
       strcmp(w, sprintf('shaolin\n')))
  if (strcmp(ELEM_TYPE, 'double'))
    prefix = 'LD_PRELOAD=/usr/lib/libfftw3.so';
  elseif (strcmp(ELEM_TYPE, 'float32'))
    prefix = 'LD_PRELOAD=/usr/lib/libfftw3f.so';
  else
    assert(0);
  end
else
  assert(0);
end
