function recovery_write_jobs(submit, only_missing)
%% =======================================================================
% Write one SLURM job file per fit (and optionally submit each with sbatch)
%
%   recovery_write_jobs            % write the job files, do not submit
%   recovery_write_jobs(true)      % write and submit
%   recovery_write_jobs(true,true) % only fits without a result file (default)
%
% Files: jobs/job_gGG_dDD_fFF.sh (each runs recovery_job(task_id))
%        logs/gGG_dDD_fFF_<jobid>.out
%% =======================================================================

if nargin < 1, submit = false; end
if nargin < 2, only_missing = true; end

%% ---------------- SLURM / MATLAB settings: EDIT THESE (or paste the header of your usual .sh files) ----------------
slurm.time      = '250:00:00';
slurm.mem       = '10G';
slurm.cpus      = 1;
slurm.partition = '';              % e.g. 'cpu'; '' = cluster default
slurm.account   = '';              % '' = none
slurm.modules   = {};              % e.g. {'matlab/R2023b'}
matlab_cmd      = '/sw/matlab2025b/bin/matlab -nodisplay -nodesktop -singleCompThread';
% extra lines to put in every job file after the #SBATCH lines (e.g. 'export MLM_LICENSE_FILE=...')
extra_lines     = {};
%% ---------------------------------------------------------------------------------------------------------------

cfg = recovery_setup();
T   = recovery_tasklist(cfg);
job_dir = fullfile(cfg.code_dir, 'jobs');
log_dir = fullfile(cfg.code_dir, 'logs');
if ~exist(job_dir,'dir'), mkdir(job_dir); end
if ~exist(log_dir,'dir'), mkdir(log_dir); end

n_written = 0;  n_skipped = 0;  n_submitted = 0;
for k = 1:size(T,1)
    g = T(k,1);  d = T(k,2);  f = T(k,3);
    tag = sprintf('g%02d_d%02d_f%02d', g, d, f);
    if only_missing && exist(fullfile(cfg.out_dir, ['rec_' tag '.mat']), 'file')
        n_skipped = n_skipped + 1;  continue
    end
    fn = fullfile(job_dir, ['job_' tag '.sh']);
    fid = fopen(fn, 'w');
    fprintf(fid, '#!/bin/sh\n');
    fprintf(fid, '#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --mail-type=NONE\n');
    fprintf(fid, '#SBATCH --job-name=rec_%s\n', tag);
    fprintf(fid, '#SBATCH --output=%s\n', fullfile(log_dir, [tag '_%j.out']));
    fprintf(fid, '#SBATCH --time=%s\n', slurm.time);
    fprintf(fid, '#SBATCH --mem=%s\n', slurm.mem);
    fprintf(fid, '#SBATCH --cpus-per-task=%d\n', slurm.cpus);
    if ~isempty(slurm.partition), fprintf(fid, '#SBATCH --partition=%s\n', slurm.partition); end
    if ~isempty(slurm.account),   fprintf(fid, '#SBATCH --account=%s\n',   slurm.account);   end
    fprintf(fid, '\n');
    for i = 1:numel(slurm.modules), fprintf(fid, 'module load %s\n', slurm.modules{i}); end
    for i = 1:numel(extra_lines),   fprintf(fid, '%s\n', extra_lines{i}); end
    fprintf(fid, 'cd "%s"\n', cfg.code_dir);
    fprintf(fid, '%s -r "addpath(''%s''); try, recovery_job(%d); catch e, disp(getReport(e)); exit(1); end; exit(0);"\n', ...
            matlab_cmd, cfg.code_dir, k);
    fclose(fid);
    fileattrib(fn, '+x');
    n_written = n_written + 1;

    if submit
        [st, out] = system(sprintf('sbatch "%s"', fn));
        if st == 0
            n_submitted = n_submitted + 1;
        else
            warning('sbatch failed for %s: %s', fn, strtrim(out));
        end
    end
end
end
