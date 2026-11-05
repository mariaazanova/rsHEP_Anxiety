function varargout = plot_spec(X,Fs,varargin)
% this function plots the psd of the input signal. The plot is interactive
% and prints the channel number by clicking on each line.

% INPUTS:
%       * X [time x channel]: the multichannel signal
%       * Fs : sampling frequency
%       is applied on the frequency limit.
%       * varargin:
%           - 'DoPlot' (bool): to plot the psd (1, default) or not (0)
%           - 'frequency_mark': the frequencies that are of specific 
%             interest. A vertical line will be ploted at those frequencies. 
%           -  'F_max': the maximum frequency of interest. 
%              If empty (default) then no limit
%           - 'freq_res': ferquency resolution (default 1Hz).
%           -  'noverlap': overlap for welch method. default=0
% OUTPUS: 
%       The outputs are assigned based on the number of outputs specified:
%       * 1st output: 
%                   - X_spec [frequency x channel]: the psd 
%       * 2nd output:
%                   - F [frequency x 1]: the frequency bins 
% ------------------------------------------------------------------------
% Usage Examples:
%          * plot_spec(X,Fs): will plot the psd
%          * plot_spec(X,Fs,'f_max',45): will plot the psd up to 45Hz
%          * plot_spec(X,Fs, 'frequency_mark', [10,20]): will plot
%          vertical dashed lines at 10Hz and 20Hz
%          * plot_spec(X,Fs, 'DoPlot', 0): does not plot, does not give
%          any output
%          * plot_spec(X,Fs, 'freq_res', 0.5): frequency resolution 0.5Hz
%          * X_spec = plot_spec(X,Fs, ... ): gives the psd in the output
%          * [X_spec,F] = plot_spec(X,Fs, ... ): gives the psd and the
%          frequency bins in the output
% ------------------------------------------------------------------------
% TODO:
%      * get the desired channels
%      * get the start and stop time
%% 
% copyright (c): Mina Jamshidi Idaji (minajamshidi91@gmail.com) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
% DAMAGE.
%% default values
freqmark = [];
DoPlot = 1;
F_max = [];
freq_res = 1;
n_overlap = 0;
%% input check
if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'doplot'
                DoPlot = varargin{j+1};
            case 'frequency_mark'
                freqmark = varargin{j+1};
            case 'f_max'
                F_max = varargin{j+1};
            case 'noverlap'
                n_overlap = varargin{j+1};
            case 'freq_res'
                freq_res = varargin{j+1};
        end
    end
end
%%
nfft = 2^ceil(log2(Fs/freq_res));
[X_spec,F] = pwelch(X,nfft,n_overlap,nfft,Fs);
if numel(F_max)
    X_spec = X_spec(F<=F_max,:);
    F = F(F<=F_max);
end
if nargout>0
    varargout{1} = X_spec;
end
if nargout>1
    varargout{2} = F;
end


if DoPlot
    hold on
    for n = 1:size(X_spec,2)
        command = [ 'disp(''x ' num2str(n) ''')' ];
        pl(n) = plot(F,10*log10(X_spec(:,n)),'linewidth',1.5, 'ButtonDownFcn', command);
    end
    ylabel('power/Hz (dB)');
    xlabel('frequency (Hz)');
    yLimit_spec = get(gca,'YLim');
    for k = 1:numel(freqmark)
        hold on,
        plot(freqmark(k)*[1,1],yLimit_spec,'--','color',rand(1,3));
    end
end

end% end function