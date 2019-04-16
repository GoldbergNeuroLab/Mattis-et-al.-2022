import os
import logging
from pyabf import ABF
import matplotlib.pyplot as plt
import numpy as np

from trace_analysis import fit_tophat, find_peaks, get_derivative
from config import ABF_LOCATION, AP_THRESHOLD
from sys import float_info

ABF_FILE_EXTENSION = '.abf'
EXPERIMENT_TYPE_CURRENT_STEPS = 'current_steps'

EXPERIMENT_TYPES = [
    EXPERIMENT_TYPE_CURRENT_STEPS
]

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class Sweep(object):
    """
    The data related to one sweep in an `ExperimentData` (i.e. in an abf file)
    """

    def __init__(
            self,
            time_steps,
            input_signal,
            output_signal,
            time_steps_units=None,
            input_signal_units=None,
            output_signal_units=None,
            sweep_name=None):
        """

        :param time_steps: The time steps of the sweep (list)
        :param input_signal: The input signal values (list)
        :param output_signal: The output signal values (list)
        :param time_steps_units: units of `time`
        :param input_signal_units: units of `input_signal`
        :param output_signal_units: units of `output_signal`
        :param metadata: Optional metadata about the sweep (dict)
        """
        # time, input signal and output signal must have the same number of
        # data points
        assert len(time_steps) == len(input_signal)
        assert len(time_steps) == len(output_signal)

        self.time_steps = time_steps
        self.input_signal = input_signal
        self.output_signal = output_signal
        self.time_steps_units = time_steps_units
        self.input_signal_units = input_signal_units
        self.output_signal_units = output_signal_units
        self.sweep_name = sweep_name
        self.analysis_cache = {}  #TODO replace this with functools.lru_cache

        logger.info('{} input units {}'.format(sweep_name, input_signal_units))
        logger.info('{} output units {}'.format(sweep_name, output_signal_units))

        # TODO cache results of analyses already completed

    def __str__(self):
        return 'Data from a single sweep, containing {} data points'.format(len(self.time_steps))

    def fit_input_tophat(
            self, verify=False, verify_file='verfication.png'):
        """
        Fit a tophat function to the input signal and cache the result

        :return: (base_level, hat_level, hat_mid, hat_width)
        """
        if 'fit_input_tophat' in self.analysis_cache:
            logger.info('Found tophat params in cache')
            return self.analysis_cache['fit_input_tophat']

        tophat_params = fit_tophat(
            self.time_steps,
            self.input_signal,
            verify=verify,
            verify_file=verify_file)
        self.analysis_cache['fit_input_tophat'] = tophat_params

        # TODO Try to spot zero height tophats, and return something sensible

        return tophat_params

    def find_output_peaks(self, threshold=0, verify=False, verify_file='verification.png'):
        """
        Find the peaks in the output and return a list of (t, V) values

        :param threshold:
        :param verify:
        :param verify_file:
        :return: list of tuples (peak time, peak value)
        """
        if 'find_output_peaks' in self.analysis_cache:
            logger.info('Found peak counts in cache')
            return self.analysis_cache['find_output_peaks']

        peaks = find_peaks(
            self.time_steps,
            self.output_signal,
            threshold=threshold,
            verify=verify,
            verify_file=verify_file)
        self.analysis_cache['find_output_peaks'] = peaks

        return peaks

    def get_output_derivative(self, verify=False):
        """
        return d/dt values of output signal. List indices correspond to
        self.time_steps

        :return: list of dV/dt values
        """
        if 'get_output_derivative' in self.analysis_cache:
            logger.info('Found output derivative in cache')
            return self.analysis_cache['get_output_derivative']

        d_dt_output = get_derivative(self.output_signal, self.time_steps)
        self.analysis_cache['get_output_derivative'] = d_dt_output

        if verify:
            plt.close('all')
            plt.figure(figsize=(16, 10))
            fig, ax1 = plt.subplots()
            ax1.plot(self.time_steps, self.output_signal)

            ax2 = ax1.twinx()
            ax2.plot(self.time_steps, d_dt_output)
            plt.show()

        return d_dt_output

    def get_output_second_derivative(self):
        """
        return d2V/dt2 values of output signal. List indices correspond to
        self.time_steps

        :return: list of d2V/dt2 values
        """
        if 'get_output_second_derivative' in self.analysis_cache:
            logger.info('Found output second derivative in cache')
            return self.analysis_cache['get_output_second_derivative']

        d2_dt2_output = get_derivative(self.get_output_derivative(), self.time_steps)
        self.analysis_cache['get_output_second_derivative'] = d2_dt2_output
        return d2_dt2_output

    def show_plot(self):
        """
        Plot input vs time and output vs time on overlapping axes
        :return: None
        """
        fig, ax1 = plt.subplots()
        ax1.plot(self.time_steps, self.input_signal, color='red')
        ax2 = ax1.twinx()
        ax2.plot(self.time_steps, self.output_signal)

        plt.show()


class ExperimentData(object):
    """The set of traces in one abf file (a colloquial, not mathematical set)"""
    def __init__(self, abf):
        """

        :param abf: The abf file, as loaded by pyabf
        :param experiment_type: The type of experiment this data is from. One of EXPERIMENT_TYPES
        """
        self.abf = abf
        self.filename = os.path.basename(abf.abfFilePath)
        self.sweep_count = abf.sweepCount
        logger.info('{} sweeps in {}'.format(self.sweep_count, self.filename))

        # Extract all the sweeps into
        self.sweeps = []
        for sweep_num in self.abf.sweepList:
            self.abf.setSweep(sweep_num, channel=1)
            input_signal = self.abf.sweepY
            input_signal_units = self.abf.sweepUnitsY

            self.abf.setSweep(sweep_num, channel=0)
            output_signal_units = self.abf.sweepUnitsY
            output_signal = self.abf.sweepY

            time_steps = self.abf.sweepX
            time_units = self.abf.sweepUnitsX
            self.sweeps.append(
                Sweep(
                    time_steps,
                    input_signal,
                    output_signal,
                    time_units,
                    input_signal_units,
                    output_signal_units,
                    sweep_name='{}_{}'.format(self.filename[:-4], sweep_num)
                )
            )

    def __str__(self):
        return('Experiment data from {} containing {} sweeps of {} data'.format(
            self.filename, self.sweep_count, self.experiment_type
        ))


class VCTestData(ExperimentData):
    """Functions to get relevant metrics for 'VC test' experiments"""
    def get_input_resistance(self):
        """
        Input resistance: calculate using change in steady state current
        in response to small hyperpolarizing voltage step

        :return:
        """
        resistances = []
        for sweep in self.sweeps:
            voltage_base, applied_voltage, voltage_mid, voltage_width = \
                sweep.fit_input_tophat()  # Voltage base is should always be ~0
            voltage_start = voltage_mid - voltage_width / 2
            logger.info('Current starts at t={}'.format(voltage_start))
            voltage_end = voltage_mid + voltage_width / 2
            start_idx = None
            end_idx = None
            for idx, t in enumerate(sweep.time_steps):
                if t > voltage_start and start_idx is None:
                    start_idx = idx
                if t > voltage_end and end_idx is None:
                    end_idx = idx

            # Measure current for the middle half of the driven part of the sweep
            logger.info('Driven slice is {} to {}'.format(start_idx, end_idx))
            measurement_slice_start = start_idx + (end_idx - start_idx) // 4
            measurement_slice_end = start_idx + 3 * (end_idx - start_idx) // 4

            mean_current_in_measurement_slice = np.mean(
                sweep.output_signal[measurement_slice_start: measurement_slice_end]
            )

            # Measure current for the middle half post drive part of the sweep
            last_idx = len(sweep.input_signal) - 1
            resting_slice_start = end_idx + (last_idx - end_idx) // 4
            resting_slice_end = end_idx + 3 * (last_idx - end_idx) // 4

            mean_current_in_resting_slice = np.mean(
                sweep.output_signal[resting_slice_start: resting_slice_end]
            )

            logger.info('Applied voltage: {} {}'.format(
                applied_voltage, sweep.input_signal_units))
            logger.info('Mean driven current: {} {}'.format(
                mean_current_in_measurement_slice, sweep.output_signal_units))
            logger.info('Resting current is: {} {}'.format(
                mean_current_in_resting_slice, sweep.input_signal_units))

            change_in_current = mean_current_in_measurement_slice - mean_current_in_resting_slice
            resistance = applied_voltage / change_in_current
            resistances.append(resistance)

        return resistances


class CurrentClampGapFreeData(ExperimentData):
    """Functions to get relevant metrics for 'current clamp gap free' experiments"""
    def get_resting_potential(self, verify=True):
        """
        Resting potential is in the output trace. Just average it. There should
        be jsut one trace

        :return:
        """
        def verification_plot():
            plt.figure(figsize=(32, 20))
            plt.close('all')
            plt.plot(self.sweeps[0].time_steps, self.sweeps[0].output_signal)
            plt.axhline(mean_voltage)
            plt.show()

        assert len(self.sweeps) == 1
        mean_voltage = np.mean(self.sweeps[0].output_signal)

        if verify:
            verification_plot()

        return mean_voltage


class CurrentStepsData(ExperimentData):
    """Functions to get relevant metrics for 'current steps' experiments"""
    def get_current_step_sizes(self):
        """
        Get a list of the step sizes of the driving current in the same order
        as self.sweeps

        :return: list of floats
        """
        step_sizes = []
        for sweep in self.sweeps:
            base_level, hat_level, hat_mid, hat_width = sweep.fit_input_tophat()
            current_step = hat_level - base_level
            step_sizes.append(current_step)

        return step_sizes

    def get_ap_counts(self):
        """
        Get a list of the number of action potentials in the output in the
        same order as self.sweeps

        :return: list of ints
        """
        ap_counts = []
        for sweep in self.sweeps:
            ap_counts.append(len(sweep.find_output_peaks()))

        return ap_counts

    def get_rheobase(self, verify=False):
        """
        Get the rheobase - the minimum voltage that elicits at least one peak
        :return:
        """
        def verification_plot():
            plt.figure(figsize=(32, 20))
            plt.close('all')
            for i, sweep in enumerate(self.sweeps):
                offset = 140 * i  # TODO Derive offset from data
                plt.plot(sweep.time_steps, sweep.output_signal + offset)
                text_height = sweep.output_signal[0] + offset
                plt.text(0, text_height, '{:.0f}'.format(peaks_at_voltage[i][1]), fontsize=8)
                if i == rheobase_sweep_num:
                    plt.text(
                        .95,
                        text_height,
                        'Rheobase from this sweep',
                        horizontalalignment='right',
                        fontsize=8
                    )

            plt.gca().get_yaxis().set_visible(False)
            plt.show()

        drive_voltages = []
        peak_count = []
        for sweep in self.sweeps:
            # Find the voltage of the driving signals for all sweeps
            tophat_params = sweep.fit_input_tophat()
            drive_voltages.append(tophat_params[1] - tophat_params[0])

            # count the peaks
            peak_count.append(len(sweep.find_output_peaks()))

        peaks_at_voltage = list(zip(peak_count, drive_voltages))

        for i, sweep in enumerate(peaks_at_voltage):
            if sweep[0] > 0:
                rheobase = sweep[1]  # return voltage of fist sweep with a peak
                rheobase_sweep_num = i
                break
        else:
            logger.info('No sweep had a peak in this experiment')
            # TODO Might be better to raise an exception
            rheobase = None  # return None if there are no APs in the experiment
            rheobase_sweep_num = None

        if verify:
            verification_plot()

        return rheobase

    def get_spike_frequency_adaptation(self):
        """
        Spike frequency adaptation: ratio of first to 10th interspike interval
        (ISI1/ISI10) and ratio of first to last interspike interval
        (ISI1/ISIn) for the first suprathreshold current injection to elicit
        sustained firing

        TODO What counts as sustained firing? First sweep with 11+ peaks?

        :return: (isi_1/isi_10, isi_1/isi_n)
        """
        peaks_by_sweep = []
        for sweep in self.sweeps:
            peaks_by_sweep.append(sweep.find_output_peaks())

        for peaks in peaks_by_sweep:
            if len(peaks) >= 11:
                isi_1 = peaks[1][0] - peaks[0][0]
                isi_10 = peaks[10][0] - peaks[9][0]
                isi_n = peaks[-1][0] - peaks[-2][0]
                return isi_1/isi_10, isi_1/isi_n
        else:
            logger.info('Data had no sweeps with sustained firing')
            return None, None

    def get_max_steady_state_firing_frequency(self, minimum_ap_amplitude=40, verify=False):
        """
         Max steady state firing frequency:
         max mean firing frequency in response to current injection with no
         failures (AP amplitude at least 40mV and overshooting 0mV)

        # TODO what should be returned. frequency. Driving voltage eliciting that frequency?
        # TODO do we have to check for "missing" peaks
        :return: frequency, inverse of timesteps units
        """
        def verification_plot():
            plt.figure(figsize=(32, 20))
            plt.close('all')
            for i, sweep in enumerate(self.sweeps):
                offset = 140 * i  # TODO Derive offset from data
                text_height = sweep.output_signal[0] + offset
                plot_color = 'b'
                for freq, idx in frequencies.items():
                    if idx == i:
                        plt.text(0, text_height, '{:.2f}Hz'.format(freq), fontsize=8)
                if i == max_frequency_idx:
                    max_freq_peaks = peaks_by_sweep[i]
                    first_ap_time = max_freq_peaks[0][0]
                    last_ap_time = max_freq_peaks[-1][0]
                    time_diff = last_ap_time - first_ap_time
                    plt.text(
                        .95,
                        text_height,
                        'Max SSFF: {} APs in {:.2f}s = {:.0f}Hz'.format(
                            len(max_freq_peaks), time_diff, max_frequency),
                        horizontalalignment='right',
                        fontsize=8
                    )
                    plt.axvline(first_ap_time)
                    plt.axvline(last_ap_time)
                    plot_color = 'r'
                plt.plot(sweep.time_steps, sweep.output_signal + offset, color=plot_color)

            plt.gca().get_yaxis().set_visible(False)
            plt.show()

        peaks_by_sweep = []
        for sweep in self.sweeps:
            # Set threshold = 0 to fulfil "overshooting 0mV" criterion
            peaks_by_sweep.append(sweep.find_output_peaks(threshold=0))

        frequencies = {}  # {frequency: sweep_num, ... }
        for i, peaks in enumerate(peaks_by_sweep):
            invalid_sweep = False
            if len(peaks) < 2:
                logger.info('Not enough peaks in sweep {} to calculate a frequency'.format(i))
                invalid_sweep = True

            for peak in peaks:
                if peak[1] - self.get_ap_threshold_1() < minimum_ap_amplitude:  # Fulfil amplitude > failed_ap_threshold
                    logger.warning('One of the peaks in sweep {} had amplitude < {} mV'.format(i, minimum_ap_amplitude))
                    invalid_sweep = True
                    break

            if not invalid_sweep:
                frequency = (len(peaks) - 1)/(peaks[-1][0] - peaks[0][0])
                frequencies[frequency] = i

        if len(frequencies) == 0:
            logger.warning('No sweep had enough peaks to calculate a frequency')
            return 0.0
        else:
            max_frequency = max(frequencies)
            max_frequency_idx = frequencies[max_frequency]
            logger.info('Max SSFF is {} from sweep {}'.format(max_frequency, max_frequency_idx))

        if verify:
            verification_plot()

        return max_frequency

    def get_max_instantaneous_firing_frequency(self):
        """
        Max instantaneous firing frequency:
        inverse of smallest interspike interval in response to current
        injection (AP amplitude at least 40mV and overshooting 0mV)


        :return:
        """
        minimum_peak_interval = float_info.max
        for sweep in self.sweeps:
            # Set threshold = 0 to fulfil "overshooting 0mV" criterion
            peaks = sweep.find_output_peaks(threshold=0)
            for idx, peak in enumerate(peaks[1:]):
                peak_interval = peak[0] - peaks[idx][0]
                logger.debug('{}, {} to {} peak interval is {}'.format(
                    sweep.sweep_name, idx, idx + 1, peak_interval
                ))
                if peak_interval < minimum_peak_interval:
                    logger.info(
                        'Found a smaller peak interval in {}, peaks {} to {}'.format(
                            sweep.sweep_name, idx, idx + 1
                        )
                    )
                    minimum_peak_interval = peak_interval

        return 1/minimum_peak_interval

    def _get_ap_threshold_1_details(self):
        """
        AP threshold #1:
        for first spike obtained at suprathreshold current injection, the
        voltage at which first derivative (dV/dt) of the AP waveform reaches
        10V/s = 10000mV/s

        :return: sweep number of threshold measurement, V or threshold, t of threshold
        """
        gradient_threshold = 10000  # V/s  TODO handle units properly

        # iterate through sweeps and peaks until we find the first peak. We will
        # return a result based on that peak.
        for sweep_num, sweep in enumerate(self.sweeps):
            for peak in sweep.find_output_peaks():
                dVdt = sweep.get_output_derivative()
                for idx, gradient in enumerate(dVdt):
                    if gradient > gradient_threshold:
                        # return the value of the voltage at the timestamp that
                        # we cross the threshold in gradient
                        logger.info('Found AP threshold 1 in {}'.format(sweep.sweep_name))
                        return sweep_num, sweep.output_signal[idx], sweep.time_steps[idx]

    def _get_ap_threshold_1_time(self):
        """


        :return:
        """
        return self._get_ap_threshold_1_details()[2]

    def get_ap_threshold_1(self):
        """

        :return:
        """
        return self._get_ap_threshold_1_details()[1]

    def get_ap_threshold_2(self):
        """
        AP threshold #2:
        for first spike obtained at suprathreshold current injection, the
        voltage at which second derivative (d2V/dt2) reaches 5% of maximal
        value

        # TODO this return value is dubious for the test data
        # TODO it occurs outside the driving voltage step.
        # TODO Would we get a reasonable result by looking at the slice in the
        # TODO driving tophat

        :return:
        """
        raise NotImplementedError
        # for sweep in self.sweeps:
        #     for peak in sweep.find_output_peaks():
        #         d2V_dt2 = sweep.get_output_second_derivative()
        #         d2V_dt2_peaks = find_peaks(sweep.time_steps, d2V_dt2)
        #         max_first_d2V_dt2_peak = d2V_dt2_peaks[0][1]
        #         for idx, d2V_dt2_value in enumerate(d2V_dt2):
        #             if d2V_dt2_value > 0.05 * max_first_d2V_dt2_peak:
        #                 logger.info('Found AP threshold 2 in {}'.format(sweep.sweep_name))
        #                 logger.info('Found AP threshold 2 at t={}'.format(sweep.time_steps[idx]))
        #                 return sweep.output_signal[idx]

    def get_ap_rise_time(self):
        """
        AP rise time:
        for first spike obtained at suprathreshold current injection, time
        from AP threshold 1 to peak

        :return:
        """
        sweep_num, ap_threshold, ap_threshold_time = self._get_ap_threshold_1_details()
        sweep = self.sweeps[sweep_num]
        ap_peak_time = sweep.find_output_peaks()[0][0]
        return ap_peak_time - ap_threshold_time

    def get_ap_amplitude(self):
        """
        AP amplitude:
        for first spike obtained at suprathreshold current injection, change
        in mV from AP threshold #1 to peak

        :return:
        """
        sweep_num, ap_threshold_voltage, ap_threshold_time = \
            self._get_ap_threshold_1_details()
        sweep = self.sweeps[sweep_num]
        ap_peak_voltage = sweep.find_output_peaks()[0][1]
        return ap_peak_voltage - ap_threshold_voltage

    def get_ap_half_width(self, verify=False):
        """
        AP half-width:
        for first spike obtained at suprathreshold current injection, width
        of the AP (in ms) at 1/2 maximal amplitude, using AP threshold #1 and
        AP amplitude

        :return:
        """
        def verification_plot():
            plt.close('all')
            plot_slice_start_idx = peak_start_idx - 3 * (peak_end_idx - peak_start_idx)
            plot_slice_end_idx = peak_end_idx + 5 * (peak_end_idx - peak_start_idx)
            plt.plot(
                sweep.time_steps[plot_slice_start_idx: plot_slice_end_idx],
                sweep.output_signal[plot_slice_start_idx: plot_slice_end_idx]
            )

            # TODO add text for threshold, peak, half, horizontal lines
            plt.axhline(threshold_amp, color='b')
            plt.axhline(sweep.output_signal[peak_idx], color='b')

            half_width_amp = 0.5 * (sweep.output_signal[peak_idx] + threshold_amp)
            plt.axhline(half_width_amp, color='r')

            plt.axvline(sweep.time_steps[peak_start_idx])
            plt.axvline(sweep.time_steps[peak_end_idx])

            plt.text(
                sweep.time_steps[plot_slice_end_idx],
                0.5 * (peak_amplitude + half_width_amp),
                'AP half width ={:.2g}'.format(ap_half_width),
                horizontalalignment='right',
                verticalalignment='top')

            plt.show()

        for sweep in self.sweeps:
            for peak in sweep.find_output_peaks():
                logger.info('Found first peak in {}'.format(sweep.sweep_name))
                threshold_amp = self.get_ap_threshold_1()
                peak_time, peak_amplitude = peak
                half_peak_voltage = 0.5 * (peak_amplitude + threshold_amp)
                peak_idx = list(sweep.time_steps).index(peak_time)
                logger.info('Peak time is {}'.format(peak_time))
                logger.info('Peak is at data point number {}'.format(peak_idx))

                voltage_at_time = dict(zip(sweep.time_steps, sweep.output_signal))
                # Iterate back through the data to find the half peak time
                for idx_diff, time_step in enumerate(sweep.time_steps[peak_idx::-1]):
                    if voltage_at_time[time_step] < half_peak_voltage:
                        logger.info('Found peak start at {}'.format(time_step))
                        peak_start = time_step
                        peak_start_idx = peak_idx - idx_diff
                        break

                # Iterate forward through the data to find the half peak time
                for idx_diff, time_step in enumerate(sweep.time_steps[peak_idx:]):
                    if voltage_at_time[time_step] < half_peak_voltage:
                        logger.info('Found peak end at {}'.format(time_step))
                        peak_end = time_step
                        peak_end_idx = peak_idx + idx_diff
                        break

                ap_half_width = peak_end - peak_start

                if verify:
                    verification_plot()

                return ap_half_width

    def get_input_resistance(self):
        """
        Input resistance #1:
        calculate using slope of the linear fit to the plot of the V-I
        relation from subthreshold current steps at/around resting potential

        # TODO do this from the last sweep before threshold?
        # TODO Average gradient for first 20ms?
        # TODO What is resting potential?
        # TODO Need to discuss this one.

        :return:
        """
        raise NotImplementedError
        # Find the sweep which elicits the first AP
        for sweep_num, sweep in enumerate(self.sweeps):
            if len(sweep.find_output_peaks()) > 0:
                first_suprathreshold_sweep = sweep_num
                break

        last_subthreshold_sweep_num = first_suprathreshold_sweep - 1

        # Find time of current step in last subthreshold peak
        last_subthreshold_sweep = self.sweeps[last_subthreshold_sweep_num]
        tophat_params = last_subthreshold_sweep.fit_input_tophat()
        current_step_time = tophat_params[2] - tophat_params[3] / 2
        drive_current = tophat_params[1]

        # Get the derivative of the output function
        dV_dt = last_subthreshold_sweep.get_output_derivative()

    def plot_v_vs_i(self, sweep_num):
        """
        Just for testing

        :param sweep_num:
        :return:
        """
        sweep = self.sweeps[sweep_num]
        fig, ax1 = plt.subplots()
        ax1.plot(sweep.output_signal, sweep.input_signal)
        plt.show()


def get_file_list(abf_location):
    """
    Figure out which file(s) to analyze based ont eh location(s) specified in the config file.
    Locations can be strings with the name of a file / folder path, or a list of strings.

    :param abf_location: The abf location as extracted from config
    :return:
    """
    if isinstance(abf_location, list):
        abf_location_list = abf_location
    else:
        abf_location_list = [abf_location]

    abf_files = []
    error = False
    for path in abf_location_list:
        if not os.path.exists(path):
            logger.error('File {} not found'.format(path))
            error = True
        if os.path.isdir(path):
            abf_files += [f for f in os.listdir(path) if f.endswith(ABF_FILE_EXTENSION)]
        elif os.path.isfile(path):
            if path.endswith(ABF_FILE_EXTENSION):
                abf_files.append(path)

    if error:
        raise ValueError('Specified location for abd files does not exist')

    logger.info('Found {} files to analyze'.format(len(abf_files)))
    return abf_files


if __name__ == '__main__':

    abf_files = get_file_list(ABF_LOCATION)
    for filename in abf_files:
        logger.info('Filename: {}'.format(filename))
        abf = ABF(filename)
        # for field in dir(abf):
        #     print('#################### {}'.format(field))
        #     exec('print(abf.{})'.format(field))

        ############################   CURRENT STEPS
        experiment = CurrentStepsData(abf)
        # for sweep in experiment.sweeps:
        #     sweep.show_plot()


        # rheobase = experiment.get_rheobase(verify=True)
        # print('Rheobase of {} is {}mV'.format(experiment.filename, rheobase))
        # exit()
        #
        # sfa = experiment.get_spike_frequency_adaptation()
        # print('SFA is {}'.format(sfa))

        max_ssff = experiment.get_max_steady_state_firing_frequency(verify=True)
        print('Max steady state firing frequency is {}'.format(max_ssff))

        # max_iff = experiment.get_max_instantaneous_firing_frequency()
        # print('Max instantaneous firing frequency is {}'.format(max_iff))
        #
        # ap_threshold_1 = experiment.get_ap_threshold_1()
        # print('AP threshold 1 is {}'.format(ap_threshold_1))
        #
        # try:
        #     ap_threshold_2 = experiment.get_ap_threshold_2()
        #     print('AP threshold 2 is {}'.format(ap_threshold_2))
        # except NotImplementedError:
        #     logger.warning("I don't know how to do that")
        #
        # ap_half_width = experiment.get_ap_half_width()
        # print('AP half width is {}'.format(ap_half_width))
        ############################   \CURRENT STEPS

        ############################   VC TEST
        # experiment = VCTestData(abf)
        # print('time units: {}, input units: {}, output units: {}'.format(
        #     experiment.sweeps[0].time_steps_units,
        #     experiment.sweeps[0].input_signal_units,
        #     experiment.sweeps[0].output_signal_units
        # ))
        # input_resistances = experiment.get_input_resistance()
        # print('Input resistances: {}'.format(input_resistances))
        # print('Input resistance is {} {}/{}'.format(
        #     np.mean(input_resistances),
        #     experiment.sweeps[0].input_signal_units,
        #     experiment.sweeps[0].output_signal_units
        # ))
        # print('mV / pA is GOhm')
        ############################   \VC TEST

        ############################   current clamp gap free
        # experiment = CurrentClampGapFreeData(abf)
        # resting_potential = experiment.get_resting_potential()
        # print('Resting potential is: {} {}'.format(
        #     resting_potential, experiment.sweeps[0].output_signal_units))
        # for sweep in experiment.sweeps:
        #     sweep.show_plot()
        #     print('Input: {}'.format(sweep.input_signal_units))
        #     print('Output: {}'.format(sweep.output_signal_units))

        ############################   /current clamp gap free