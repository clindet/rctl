var assert = require('assert')
var expect = require('chai').expect;
var should = require('chai').should();
var rbashful = require('./rbashful_exec');

describe('Test of command line tool rbashful', function() {
    it("should contain No commands were ran", function(){
        text = rbashful();
        expectStr = "No commands were ran";
        assert.equal(text.indexOf(expectStr) > 0, true);
    })
})
