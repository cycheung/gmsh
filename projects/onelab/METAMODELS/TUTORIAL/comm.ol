LOGFILES.radioButton(0);
SHOWSOCKETMESSAGES.radioButton(0);

MESSAGE.string(Hello World!);
OL.show(MESSAGE);

Interfaced.register(interfaced);
Interfaced.run(OL.get(MESSAGE));
Interfaced.active(1);

Encapsulated.register(encapsulated);
Encapsulated.run(OL.get(MESSAGE));
Encapsulated.active(1);

Native.register(native);
Native.run();

