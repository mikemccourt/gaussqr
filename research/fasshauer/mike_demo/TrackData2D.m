%-------------------------------------------------------------------------%
%
% File: TrackData2D(ntracks,ntrackpoints)
%
% Goal: script that generates track data in 2D
%
% Inputs:   ntracks:      number of tracks
%           ntrackpoints: number of points (dsites) on each track
%
% Outputs:  N:            number of track data (dsites)
%           dsites:       NX2 matrix representing a set of N data sites
%           ylines:       vector with y coordinates of lines on tracks
%
%-------------------------------------------------------------------------%
function [N dsites ylines] = TrackData2D(ntracks,ntrackpoints)

% fix seed to generate repeatable random points
rand('seed',10);

% number of points (dsites)
N = ntracks*ntrackpoints;

% distance of lines between tracks
delta = 1/ntracks;

const = 1/ntrackpoints;

% xdsites on tracks
for j = 0:(ntracks-1)
    for i = 1:ntrackpoints
        a(i) = const*(i-1); b(i) = const*i;
        xdsites(j*ntrackpoints+i) = a(i) + rand*abs(b(i)-a(i));
    end
end

% lines on tracks
 for i = 1:ntracks
     ylines(i) = (2*i-1)*delta/2;
 end

% ydsites on tracks
for j = 0:(ntracks-1)
    for h = 1:ntrackpoints
        ydsites(h+j*ntrackpoints) = ylines(j+1)+ ...
            rand*((delta/2)-2/(ntrackpoints^2*delta)) - ... 
            (delta/4) + 1/(ntrackpoints^2*delta);
    end
end

dsites = [xdsites' ydsites'];

% % plot of points
% [xdsites' ydsites'];
% plot(xdsites',ydsites','.')
% axis square