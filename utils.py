import numpy as np
from scipy import signal
from ecg_detectors.ecgdetectors import searchBack


def hamilton_detector(unfiltered_ecg, fs, filtered_ecg):

    diff = abs(np.diff(filtered_ecg))

    b = np.ones(int(0.08 * fs))
    b = b / int(0.08 * fs)
    a = [1]

    ma = signal.lfilter(b, a, diff)

    ma[0:len(b) * 2] = 0

    peaks, _ = signal.find_peaks(ma, distance=(0.25 * fs))

    n_pks = []
    n_pks_ave = 0.0
    s_pks = []
    s_pks_ave = 0.0
    QRS = []
    RR = []
    RR_ave = 0.0

    th = 0.0

    i = 0
    idx = []
    for peak in peaks:

        if ma[peak] > th:
            QRS.append(peak)
            idx.append(i)
            s_pks.append(ma[peak])
            if len(n_pks) > 8:
                s_pks.pop(0)
            s_pks_ave = np.mean(s_pks)

            if RR_ave != 0.0:
                if QRS[-1] - QRS[-2] > 1.5 * RR_ave:
                    missed_peaks = peaks[idx[-2] + 1:idx[-1]]
                    for missed_peak in missed_peaks:
                        if missed_peak - peaks[idx[-2]] > int(0.360 * fs) and ma[missed_peak] > 0.5 * th:
                            QRS.append(missed_peak)
                            QRS.sort()
                            break

            if len(QRS) > 2:
                RR.append(QRS[-1] - QRS[-2])
                if len(RR) > 8:
                    RR.pop(0)
                RR_ave = int(np.mean(RR))

        else:
            n_pks.append(ma[peak])
            if len(n_pks) > 8:
                n_pks.pop(0)
            n_pks_ave = np.mean(n_pks)

        th = n_pks_ave + 0.45 * (s_pks_ave - n_pks_ave)

        i += 1

    QRS.pop(0)

    window = int(0.1 * fs)

    r_peaks = searchBack(QRS, unfiltered_ecg, window)

    return r_peaks


def christov_detector(fs, unfiltered_ecg):
    total_taps = 0

    b = np.ones(int(0.02*fs))
    b = b/int(0.02*fs)
    total_taps += len(b)
    a = [1]

    MA1 = signal.lfilter(b, a, unfiltered_ecg)

    b = np.ones(int(0.028*fs))
    b = b/int(0.028*fs)
    total_taps += len(b)
    a = [1]

    MA2 = signal.lfilter(b, a, MA1)

    Y = []
    for i in range(1, len(MA2)-1):

        diff = abs(MA2[i+1]-MA2[i-1])

        Y.append(diff)

    b = np.ones(int(0.040*fs))
    b = b/int(0.040*fs)
    total_taps += len(b)
    a = [1]

    MA3 = signal.lfilter(b, a, Y)

    MA3[0:total_taps] = 0

    ms50 = int(0.05*fs)
    ms200 = int(0.2*fs)
    ms1200 = int(1.2*fs)
    ms350 = int(0.35*fs)

    M = 0
    newM5 = 0
    M_list = []
    MM = []
    M_slope = np.linspace(1.0, 0.6, ms1200-ms200)
    F = 0
    F_list = []
    R = 0
    RR = []
    Rm = 0
    R_list = []

    MFR = 0
    MFR_list = []

    QRS = []

    for i in range(len(MA3)):

        # M
        if i < 5*fs:
            M = 0.6*np.max(MA3[:i+1])
            MM.append(M)
            if len(MM)>5:
                MM.pop(0)

        elif QRS and i < QRS[-1]+ms200:
            newM5 = 0.6*np.max(MA3[QRS[-1]:i])
            if newM5>1.5*MM[-1]:
                newM5 = 1.1*MM[-1]

        elif QRS and i == QRS[-1]+ms200:
            if newM5==0:
                newM5 = MM[-1]
            MM.append(newM5)
            if len(MM)>5:
                MM.pop(0)
            M = np.mean(MM)

        elif QRS and i > QRS[-1]+ms200 and i < QRS[-1]+ms1200:

            M = np.mean(MM)*M_slope[i-(QRS[-1]+ms200)]

        elif QRS and i > QRS[-1]+ms1200:
            M = 0.6*np.mean(MM)

        # F
        if i > ms350:
            F_section = MA3[i-ms350:i]
            max_latest = np.max(F_section[-ms50:])
            max_earliest = np.max(F_section[:ms50])
            F = F + ((max_latest-max_earliest)/150.0)

        # R
        if QRS and i < QRS[-1]+int((2.0/3.0*Rm)):

            R = 0

        elif QRS and i > QRS[-1]+int((2.0/3.0*Rm)) and i < QRS[-1]+Rm:

            dec = (M-np.mean(MM))/1.4
            R = 0 + dec


        MFR = M+F+R
        M_list.append(M)
        F_list.append(F)
        R_list.append(R)
        MFR_list.append(MFR)

        if not QRS and MA3[i]>MFR:
            QRS.append(i)

        elif QRS and i > QRS[-1]+ms200 and MA3[i]>MFR:
            QRS.append(i)
            if len(QRS)>2:
                RR.append(QRS[-1]-QRS[-2])
                if len(RR)>5:
                    RR.pop(0)
                Rm = int(np.mean(RR))

    QRS.pop(0)
    r_peaks = []
    search_samples = int(0.05*fs)

    for i in QRS:
        if i < search_samples:
            section = unfiltered_ecg[0:search_samples]
            r_peaks.append(np.argmax(section))
        elif i+search_samples>len(unfiltered_ecg):
            section = unfiltered_ecg[i:]
            r_peaks.append(np.argmax(section)+i)
        else:
            section = unfiltered_ecg[i-search_samples:i+search_samples]
            r_peaks.append(np.argmax(section)+(i-search_samples))

    return r_peaks


def findpeaks(data, spacing=1, limit=None):
    """
    Janko Slavic peak detection algorithm and implementation.
        https://github.com/jankoslavic/py-tools/tree/master/findpeaks

    Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """
    len = data.size
    x = np.zeros(len + 2 * spacing)
    x[:spacing] = data[0] - 1.e-6
    x[-spacing:] = data[-1] - 1.e-6
    x[spacing:spacing + len] = data
    peak_candidate = np.zeros(len)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start: start + len]  # before
        start = spacing
        h_c = x[start: start + len]  # central
        start = spacing + s + 1
        h_a = x[start: start + len]  # after
        peak_candidate = np.logical_and(peak_candidate, np.logical_and(h_c > h_b, h_c > h_a))

    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
    if limit is not None:
        ind = ind[data[ind] > limit]
    return ind


def engzee_detector(filtered_ecg, fs, unfiltered_ecg):

    diff = np.zeros(len(filtered_ecg))
    for i in range(4, len(diff)):
        diff[i] = filtered_ecg[i]-filtered_ecg[i-4]

    ci = [1,4,6,4,1]
    low_pass = signal.lfilter(ci, 1, diff)

    low_pass[:int(0.2*fs)] = 0

    ms200 = int(0.2*fs)
    ms1200 = int(1.2*fs)
    ms160 = int(0.16*fs)
    neg_threshold = int(0.01*fs)

    M = 0
    M_list = []
    neg_m = []
    MM = []
    M_slope = np.linspace(1.0, 0.6, ms1200-ms200)

    QRS = []
    r_peaks = []

    counter = 0

    thi_list = []
    thi = False
    thf_list = []
    thf = False

    for i in range(len(low_pass)):

        # M
        if i < 5*fs:
            M = 0.6*np.max(low_pass[:i+1])
            MM.append(M)
            if len(MM)>5:
                MM.pop(0)

        elif QRS and i < QRS[-1]+ms200:

            newM5 = 0.6*np.max(low_pass[QRS[-1]:i])

            if newM5>1.5*MM[-1]:
                newM5 = 1.1*MM[-1]

        elif QRS and i == QRS[-1]+ms200:
            MM.append(newM5)
            if len(MM)>5:
                MM.pop(0)
            M = np.mean(MM)

        elif QRS and i > QRS[-1]+ms200 and i < QRS[-1]+ms1200:

            M = np.mean(MM)*M_slope[i-(QRS[-1]+ms200)]

        elif QRS and i > QRS[-1]+ms1200:
            M = 0.6*np.mean(MM)

        M_list.append(M)
        neg_m.append(-M)


        if not QRS and low_pass[i]>M:
            QRS.append(i)
            thi_list.append(i)
            thi = True

        elif QRS and i > QRS[-1]+ms200 and low_pass[i]>M:
            QRS.append(i)
            thi_list.append(i)
            thi = True

        if thi and i<thi_list[-1]+ms160:
            if low_pass[i]<-M and low_pass[i-1]>-M:
                #thf_list.append(i)
                thf = True

            if thf and low_pass[i]<-M:
                thf_list.append(i)
                counter += 1

            elif low_pass[i]>-M and thf:
                counter = 0
                thi = False
                thf = False

        elif thi and i>thi_list[-1]+ms160:
                counter = 0
                thi = False
                thf = False

        if counter>neg_threshold:
            unfiltered_section = unfiltered_ecg[thi_list[-1]-int(0.01*fs):i]
            r_peaks.append(np.argmax(unfiltered_section)+thi_list[-1]-int(0.01*fs))
            counter = 0
            thi = False
            thf = False

    return r_peaks
