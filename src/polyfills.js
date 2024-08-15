import { Buffer } from 'buffer';
import stream from 'stream-browserify';
import process from 'process/browser';
import assert from 'assert';

global.Buffer = Buffer;
global.process = process;
global.stream = stream;
global.assert = assert;
