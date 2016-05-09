///////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2016 F4EXB                                                      //
// written by Edouard Griffiths                                                  //
//                                                                               //
// This program is free software; you can redistribute it and/or modify          //
// it under the terms of the GNU General Public License as published by          //
// the Free Software Foundation as version 3 of the License, or                  //
//                                                                               //
// This program is distributed in the hope that it will be useful,               //
// but WITHOUT ANY WARRANTY; without even the implied warranty of                //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                  //
// GNU General Public License V3 for more details.                               //
//                                                                               //
// You should have received a copy of the GNU General Public License             //
// along with this program. If not, see <http://www.gnu.org/licenses/>.          //
///////////////////////////////////////////////////////////////////////////////////

#include <unistd.h>

#include "dsp/dvserialworker.h"
#include "audio/audiofifo.h"

MESSAGE_CLASS_DEFINITION(DVSerialWorker::MsgMbeDecode, Message)
MESSAGE_CLASS_DEFINITION(DVSerialWorker::MsgTest, Message)

DVSerialWorker::DVSerialWorker() :
    m_running(false),
    m_currentGainIn(0),
    m_currentGainOut(0),
    m_upsamplerLastValue(0),
    m_phase(0)
{
    m_audioBuffer.resize(48000);
    m_audioBufferFill = 0;
}

DVSerialWorker::~DVSerialWorker()
{
}

bool DVSerialWorker::open(const std::string& serialDevice)
{
    return m_dvController.open(serialDevice);
}

void DVSerialWorker::close()
{
    m_dvController.close();
}

void DVSerialWorker::process()
{
    m_running  = true;
    qDebug("DVSerialWorker::process: started");

    while (m_running)
    {
        sleep(1);
    }

    qDebug("DVSerialWorker::process: stopped");
    emit finished();
}

void DVSerialWorker::stop()
{
    m_running = false;
}

void DVSerialWorker::handleInputMessages()
{
    Message* message;

    while ((message = m_inputMessageQueue.pop()) != 0)
    {
        if (MsgMbeDecode::match(*message))
        {
            MsgMbeDecode *decodeMsg = (MsgMbeDecode *) message;
            int dBVolume = (decodeMsg->getVolumeIndex() - 50) / 2;

            if (m_dvController.decode(m_dvAudioSamples, decodeMsg->getMbeFrame(), decodeMsg->getMbeRate(), dBVolume))
            {
                upsample6(m_dvAudioSamples, SerialDV::MBE_AUDIO_BLOCK_SIZE, decodeMsg->getAudioFifo());
//                upsample6(m_dvAudioSamples, m_audioSamples, SerialDV::MBE_AUDIO_BLOCK_SIZE);
//                decodeMsg->getAudioFifo()->write((const quint8 *) m_audioSamples, SerialDV::MBE_AUDIO_BLOCK_SIZE * 6, 10);
            }
            else
            {
                qDebug("DVSerialWorker::handleInputMessages: MsgMbeDecode: decode failed");
            }
        }

        delete message;
    }
}

void DVSerialWorker::upsample6(short *in, int nbSamplesIn, AudioFifo *audioFifo)
{
    for (int i = 0; i < nbSamplesIn; i++)
    {
        int cur = (int) in[i];
        int prev = (int) m_upsamplerLastValue;
        qint16 upsample;

        for (int j = 1; j < 7; j++)
        {
            upsample = (qint16) ((cur*j + prev*(6-j)) / 6);
            m_audioBuffer[m_audioBufferFill].l = upsample;
            m_audioBuffer[m_audioBufferFill].r = upsample;
            ++m_audioBufferFill;

            if (m_audioBufferFill >= m_audioBuffer.size())
            {
                uint res = audioFifo->write((const quint8*)&m_audioBuffer[0], m_audioBufferFill, 10);

                if (res != m_audioBufferFill)
                {
                    qDebug("DVSerialWorker::upsample6: %u/%u audio samples written", res, m_audioBufferFill);
                }

                m_audioBufferFill = 0;
            }
        }

        m_upsamplerLastValue = in[i];
    }
}

void DVSerialWorker::upsample6(short *in, short *out, int nbSamplesIn)
{
    for (int i = 0; i < nbSamplesIn; i++)
    {
        int cur = (int) in[i];
        int prev = (int) m_upsamplerLastValue;
        short up;

//        for (int j = 0; j < 6; j++)
//        {
//            up = 32768.0f * cos(m_phase);
//            *out = up;
//            out ++;
//            *out = up;
//            out ++;
//            m_phase += M_PI / 6.0;
//        }
//
//        if ((i % 2) == 1)
//        {
//            m_phase = 0.0f;
//        }

        up = (cur*1 + prev*5) / 6;
        *out = up;
        out++;
        *out = up;
        out++;

        up = (cur*2 + prev*4) / 6;
        *out = up;
        out++;
        *out = up;
        out++;

        up = (cur*3 + prev*3) / 6;
        *out = up;
        out++;
        *out = up;
        out++;

        up = (cur*4 + prev*2) / 6;
        *out = up;
        out++;
        *out = up;
        out++;

        up = (cur*5 + prev*1) / 6;
        *out = up;
        out++;
        *out = up;
        out++;

        up = in[i];
        *out = up;
        out++;
        *out = up;
        out++;

        m_upsamplerLastValue = in[i];
    }
}
