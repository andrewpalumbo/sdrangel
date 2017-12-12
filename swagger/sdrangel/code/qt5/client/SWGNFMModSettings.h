/**
 * SDRangel
 * This is the web API of SDRangel SDR software. SDRangel is an Open Source Qt5/OpenGL 3.0+ GUI and server Software Defined Radio and signal analyzer in software. It supports Airspy, BladeRF, HackRF, LimeSDR, PlutoSDR, RTL-SDR, SDRplay RSP1 and FunCube
 *
 * OpenAPI spec version: 4.0.0
 * Contact: f4exb06@gmail.com
 *
 * NOTE: This class is auto generated by the swagger code generator program.
 * https://github.com/swagger-api/swagger-codegen.git
 * Do not edit the class manually.
 */

/*
 * SWGNFMModSettings.h
 * 
 * 
 */

#ifndef SWGNFMModSettings_H_
#define SWGNFMModSettings_H_

#include <QJsonObject>


#include "SWGCWKeyerSettings.h"
#include <QString>

#include "SWGObject.h"


namespace SWGSDRangel {

class SWGNFMModSettings: public SWGObject {
public:
    SWGNFMModSettings();
    SWGNFMModSettings(QString* json);
    virtual ~SWGNFMModSettings();
    void init();
    void cleanup();

    QString asJson ();
    QJsonObject* asJsonObject();
    void fromJsonObject(QJsonObject &json);
    SWGNFMModSettings* fromJson(QString &jsonString);

    qint32 getBasebandSampleRate();
    void setBasebandSampleRate(qint32 baseband_sample_rate);

    qint32 getOutputSampleRate();
    void setOutputSampleRate(qint32 output_sample_rate);

    qint64 getInputFrequencyOffset();
    void setInputFrequencyOffset(qint64 input_frequency_offset);

    float getRfBandwidth();
    void setRfBandwidth(float rf_bandwidth);

    float getAfBandwidth();
    void setAfBandwidth(float af_bandwidth);

    float getFmDeviation();
    void setFmDeviation(float fm_deviation);

    float getToneFrequency();
    void setToneFrequency(float tone_frequency);

    float getVolumeFactor();
    void setVolumeFactor(float volume_factor);

    qint32 getAudioSampleRate();
    void setAudioSampleRate(qint32 audio_sample_rate);

    qint32 getChannelMute();
    void setChannelMute(qint32 channel_mute);

    qint32 getPlayLoop();
    void setPlayLoop(qint32 play_loop);

    qint32 getCtcssOn();
    void setCtcssOn(qint32 ctcss_on);

    qint32 getCtcssIndex();
    void setCtcssIndex(qint32 ctcss_index);

    qint32 getRgbColor();
    void setRgbColor(qint32 rgb_color);

    QString* getTitle();
    void setTitle(QString* title);

    qint32 getModAfInput();
    void setModAfInput(qint32 mod_af_input);

    SWGCWKeyerSettings* getCwKeyer();
    void setCwKeyer(SWGCWKeyerSettings* cw_keyer);


private:
    qint32 baseband_sample_rate;
    qint32 output_sample_rate;
    qint64 input_frequency_offset;
    float rf_bandwidth;
    float af_bandwidth;
    float fm_deviation;
    float tone_frequency;
    float volume_factor;
    qint32 audio_sample_rate;
    qint32 channel_mute;
    qint32 play_loop;
    qint32 ctcss_on;
    qint32 ctcss_index;
    qint32 rgb_color;
    QString* title;
    qint32 mod_af_input;
    SWGCWKeyerSettings* cw_keyer;
};

}

#endif /* SWGNFMModSettings_H_ */
