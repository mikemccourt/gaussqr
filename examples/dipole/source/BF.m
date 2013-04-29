function BF = BF( evalpnts,srcpnts,dipmom )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

M = size(evalpnts,1);
BF = zeros(M,3);
for i = 1:M
        R = evalpnts(i,:) - srcpnts;
        R_norm = R(:,1)*R(:,1) + R(:,2)*R(:,2) + R(:,3)*R(:,3);
        R_norm = R_norm.^(3/2);
        d = [dipmom(:,2)*R(:,3) - dipmom(:,3)*R(:,2), ...
             dipmom(:,3)*R(:,1) - dipmom(:,1)*R(:,3), ...
             dipmom(:,1)*R(:,2) - dipmom(:,2)*R(:,1)];
        d = bsxfun(@rdivide,d,R_norm);
        BF(i,:) = sum(d,1);
end
BF = 1e-07*BF;
end

% function BF = BF( evalpnts,srcpnts,dipmom )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% 
% M = size(evalpnts,1); N = size(srcpnts,1);
% BF = zeros(M,3);
% for i = 1:M
%     for j = 1:N
%         R = evalpnts(i,:) - srcpnts(j,:);
%         R_norm = R(1)*R(1) + R(2)*R(2) + R(3)*R(3);
%         R_norm = R_norm^(3/2);
%         d = [dipmom(j,2)*R(3) - dipmom(j,3)*R(2), ...
%              dipmom(j,3)*R(1) - dipmom(j,1)*R(3), ...
%              dipmom(j,1)*R(2) - dipmom(j,2)*R(1)];
%         BF(i,:) = BF(i,:) + d/R_norm;
%     end
% end
% BF = 1e-07*BF;
% end