%% Author: Muhammad Asad Ullah
%% University of Oulu

SF=[7,8,9,10,11,12];
npreamble=8;
CR = 1;
Overhead = 5;

PL = Application_payload + Overhead;

for i=1:1:length(SF)
    Ts=(2.^SF(i))/BW;
    Tpreamble=(npreamble+4.24)*Ts;

    if(LDRO==0)
        C= ceil((8.*PL-4.*SF(i)+28+16.*1)/(4*SF(i)));
    else
       % LRDO enabled: reduce bits per symbol by two as SF(i)-2
       C = ceil((8.*PL-4.*SF(i)+28+16.*1)/(4*(SF(i)-2))); 
    end

    N_payload_symbols=8+max(C.*(4+CR),0);

    Time_payload=N_payload_symbols*Ts;
    ToA(i,:)=round((Tpreamble+Time_payload).*1e3); % milliseconds
end
