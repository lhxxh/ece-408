function output_cfg_file_1D(cfg)

fid = fopen(cfg.cfg_filename, 'w');
assert(fid ~= -1);

fprintf(fid, 'dir = %s\n', cfg.dir_name);

fprintf(fid, 'x_hat_filename = %s\n', cfg.x_hat_filename);
if (cfg.save_trace)
  fprintf(fid, 'save_trace = true\n');
  fprintf(fid, 'trace_filename = %s\n', cfg.trace_filename);
else
  fprintf(fid, 'save_trace = false\n');
end

fprintf(fid, 'suffix = %s\n', cfg.suffix);

fprintf(fid, 'suffix_filename = %s\n', cfg.suffix_filename);

fprintf(fid, 'rank = %d\n', cfg.rank);
fprintf(fid, 'n = {');
for i = 1:cfg.rank-1
  fprintf(fid, '%d, ', cfg.n(i));
end
fprintf(fid, '%d}\n', cfg.n(cfg.rank));

fprintf(fid, 'N = %d\n', cfg.N);
fprintf(fid, 'T = %d\n', cfg.T);
fprintf(fid, 'L = %d\n', cfg.L);
fprintf(fid, 'M_block = %d\n', cfg.M_block);

if (cfg.quiet_mode)
  fprintf(fid, 'quiet_mode = true\n');
else
  fprintf(fid, 'quiet_mode = false\n');
end

fprintf(fid, 'lambda = %f\n', cfg.lambda);

if (cfg.poisson_noise)
  fprintf(fid, 'poisson_noise = true\n');
  fprintf(fid, 'poisson_eps = %f\n', cfg.poisson_eps);
else
  fprintf(fid, 'poisson_noise = false\n');
end

fprintf(fid, 'F_equal_I = true\n');

if (cfg.randn_debug)
  fprintf(fid, 'randn_debug = true\n');
else
  fprintf(fid, 'randn_debug = false\n');
end

if (cfg.lenkf_debug)
  fprintf(fid, 'lenkf_debug = true');
else
  fprintf(fid, 'lenkf_debug = false');
end

fclose(fid);


fid = fopen([cfg.dir_name, '/', cfg.suffix_filename], 'w');
assert(fid ~= -1);

for i=1:cfg.T
    fprintf(fid, '%s_%d\n', cfg.suffix, i-1);
end

fclose(fid);