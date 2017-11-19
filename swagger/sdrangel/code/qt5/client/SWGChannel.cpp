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


#include "SWGChannel.h"

#include "SWGHelpers.h"

#include <QJsonDocument>
#include <QJsonArray>
#include <QObject>
#include <QDebug>

namespace Swagger {

SWGChannel::SWGChannel(QString* json) {
    init();
    this->fromJson(*json);
}

SWGChannel::SWGChannel() {
    init();
}

SWGChannel::~SWGChannel() {
    this->cleanup();
}

void
SWGChannel::init() {
    index = 0;
    id = new QString("");
    title = new QString("");
    delta_frequency = 0;
}

void
SWGChannel::cleanup() {
    

    if(id != nullptr) {
        delete id;
    }

    if(title != nullptr) {
        delete title;
    }

}

SWGChannel*
SWGChannel::fromJson(QString &json) {
    QByteArray array (json.toStdString().c_str());
    QJsonDocument doc = QJsonDocument::fromJson(array);
    QJsonObject jsonObject = doc.object();
    this->fromJsonObject(jsonObject);
    return this;
}

void
SWGChannel::fromJsonObject(QJsonObject &pJson) {
    ::Swagger::setValue(&index, pJson["index"], "qint32", "");
    ::Swagger::setValue(&id, pJson["id"], "QString", "QString");
    ::Swagger::setValue(&title, pJson["title"], "QString", "QString");
    ::Swagger::setValue(&delta_frequency, pJson["deltaFrequency"], "qint32", "");
}

QString
SWGChannel::asJson ()
{
    QJsonObject* obj = this->asJsonObject();
    
    QJsonDocument doc(*obj);
    QByteArray bytes = doc.toJson();
    return QString(bytes);
}

QJsonObject*
SWGChannel::asJsonObject() {
    QJsonObject* obj = new QJsonObject();
    
    obj->insert("index", QJsonValue(index));

    toJsonValue(QString("id"), id, obj, QString("QString"));

    toJsonValue(QString("title"), title, obj, QString("QString"));

    obj->insert("deltaFrequency", QJsonValue(delta_frequency));

    return obj;
}

qint32
SWGChannel::getIndex() {
    return index;
}
void
SWGChannel::setIndex(qint32 index) {
    this->index = index;
}

QString*
SWGChannel::getId() {
    return id;
}
void
SWGChannel::setId(QString* id) {
    this->id = id;
}

QString*
SWGChannel::getTitle() {
    return title;
}
void
SWGChannel::setTitle(QString* title) {
    this->title = title;
}

qint32
SWGChannel::getDeltaFrequency() {
    return delta_frequency;
}
void
SWGChannel::setDeltaFrequency(qint32 delta_frequency) {
    this->delta_frequency = delta_frequency;
}


}

